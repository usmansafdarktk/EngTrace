import math
import random
from decimal import Decimal, ROUND_HALF_UP


def _hu(x, places):
    """Round half-up to `places` dp, resolving the tie in DECIMAL.

    `round()` resolves a half-way tie on the binary value, so
    `round(0.01185 * 1000, 1)` gives 11.8 where a reader doing decimal
    arithmetic gets 11.9 (spec P2 as amended, DECISIONS D-012). Accepts a
    Decimal so an exact decimal product can be quantised without a detour
    through a float. Matches the industrial branch's `_hu` convention.
    """
    q = Decimal(1).scaleb(-places)
    d = x if isinstance(x, Decimal) else Decimal(repr(x))
    v = d.quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


def _is_display_tie(x, places, rel_band=1e-12):
    """Is `x` at, or within a hair of, a half-way tie at `places` dp?

    A tie is the one case where NO rounding convention is defensible. A reader
    doing decimal arithmetic and applying half-up reads 0.02325 m as 23.3 mm; a
    reader using binary floats and `round()` reads it as 23.2 mm; and the
    printed line closes for exactly one of them whichever the template picks.
    Breaking the tie in decimal (P2 as amended) does not remove the ambiguity,
    it only moves it to the other reader - which is why such instances are
    RESAMPLED rather than resolved (D-016).

    Testing for an EXACT tie is not enough. Two independent evaluations of the
    same exact quantity - this template's, in kN and kN/m^2, and a solver's, in
    N and Pa - differ by a few ulps, so a rational tie such as 29.25 mm lands
    as 29.249999999999996 on one side and 29.250000000000004 on the other.
    Neither is exactly a tie, and the two then round in opposite directions:
    that is how eight non-closing instances per 60,000 seeds survived the first
    version of this guard (Phase 1 review A, finding F-2).

    So a narrow BAND is quarantined rather than a point. The band is a few
    thousand ulps wide where the tie lattice is 10^-places apart, so it removes
    nothing that is not genuinely ambiguous.
    """
    scaled = abs(x) * 10.0 ** places
    band = max(scaled * rel_band, 1e-9)
    return abs((scaled - math.floor(scaled)) - 0.5) <= band


# Template 15 (Intermediate) — Area A2: Influence Lines
def template_influence_line_max_reaction():
    """
    Maximum Support Reaction Under a Moving Two-Axle Load

    Scenario:
        A two-axle vehicle crosses a simply supported beam. The influence
        line for the reaction at A is linear (1 at A, 0 at B), so the
        maximum reaction occurs with the heavier axle directly over A and
        the lighter axle a fixed axle-spacing s inside the span:

            y(x) = (L - x) / L
            R_A,max = P_heavy * 1 + P_light * (L - s)/L

        Which axle is heavier is sampled, so the governing arrangement
        (vehicle facing left vs right) is a positioning decision.

    Difficulty: Intermediate
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 6
        (Influence Lines for Statically Determinate Structures) —
        reaction influence line and maximum-effect positioning of a
        wheel-load series.
    Physical bounds: span L in [8.0, 16.0] m; axle loads in [40, 160] kN
        differing by at least 10 kN; axle spacing s in [2.0, 4.5] m and
        always < L/3; R_max strictly between the heavy axle and the axle
        sum.

    Trace integrity (Layer 0, 2026-09-23):
        The maximum reaction R_A,max = P_heavy*1 + P_light*y (an integer
        times a 4-dp ordinate) is exact at 4 dp and was printed at 2 dp,
        where it sat on a half-way tie on ~2.6% of draws (77 + 50 * 0.8701
        = 120.5050). R_A,max is the answer, so on the owner's decision of
        2026-09-23 (D-044) its display is lengthened to 4 dp rather than
        the draw resampled: it is bound half-up in decimal at 4 dp and
        printed at 4 dp in Step 4 and in the answer. The light-axle
        ordinate y = (L - s)/L is a quotient by an arbitrary 1-dp span,
        exact at no fixed display, so a draw whose y lands on a 4-dp tie is
        resampled rather than rounded either way (D-016/D-037).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (heavier-axle identity is the positioning branch).
    # Bounded redraw: a draw whose light-axle ordinate lands on a 4-dp
    # half-way tie has no defensible gold answer and is rejected (D-016).
    for _attempt in range(200):
        L = round(random.uniform(8.0, 16.0), 1)
        P1 = random.randint(40, 160)                  # front axle
        P2 = random.choice([p for p in range(40, 161)
                            if abs(p - P1) >= 10])    # rear axle
        s = round(random.uniform(2.0, min(4.5, L / 3 - 0.1)), 1)

        P_heavy, P_light = (P1, P2) if P1 > P2 else (P2, P1)
        heavy_name = "front" if P1 > P2 else "rear"

        # 2. Core computation — round-then-recompute at every step.
        y_light_exact = (L - s) / L
        if _is_display_tie(y_light_exact, 4):
            continue                    # no defensible gold answer; redraw
        y_light = round(y_light_exact, 4)
        # An integer times a 4-dp ordinate is exact at 4 dp; bound half-up
        # in decimal and printed at 4 dp (D-044).
        R_max = _hu(
            Decimal(P_heavy) + Decimal(P_light) * Decimal(repr(y_light)), 4)
        break
    else:
        raise RuntimeError(
            "influence_line_max_reaction: no closing sample in 200 draws")

    assert s < L / 3, f"axle spacing too large: {s} vs L={L}"
    assert P_heavy < R_max < P1 + P2, f"R_max out of range: {R_max}"

    # 3. Serialize.
    question = (
        f"A simply supported beam AB (pin at A, roller at B) has a span "
        f"of L = {L:.1f} m. A two-axle vehicle crosses the beam: the "
        f"front axle carries {P1} kN and the rear axle carries {P2} kN, "
        f"with a fixed axle spacing of s = {s:.1f} m. The vehicle may "
        f"travel across the beam in either direction, and both axles "
        f"remain on the span. Using the influence line for the vertical "
        f"reaction at A, determine the maximum value of that reaction "
        f"in kN."
    )

    solution = (
        f"**Given:**\n"
        f"Span (L): {L:.1f} m\n"
        f"Axle loads: front {P1} kN, rear {P2} kN; spacing s = {s:.1f} m\n\n"
        f"**Step 1:** Construct the influence line for the reaction at "
        f"A.\n"
        f"For a unit load at distance x from A, R_A = (L - x)/L: the "
        f"ordinate is 1 at A (x = 0) and decreases linearly to 0 at B "
        f"(x = L).\n\n"
        f"**Step 2:** Position the axles for the maximum effect.\n"
        f"Since the ordinates decrease away from A, the heavier axle "
        f"(the {heavy_name} axle, {P_heavy} kN) is placed directly over "
        f"A, and the lighter axle ({P_light} kN) then sits at "
        f"x = s = {s:.1f} m.\n\n"
        f"**Step 3:** Read the influence-line ordinates under each "
        f"axle.\n"
        f"Under the heavy axle: y = 1.0000 (at A)\n"
        f"Under the light axle: y = (L - s)/L = ({L:.1f} - {s:.1f}) / "
        f"{L:.1f} = {y_light:.4f}\n\n"
        f"**Step 4:** Sum load times ordinate.\n"
        f"R_A,max = {P_heavy} * 1.0000 + {P_light} * {y_light:.4f} "
        f"= {R_max:.4f} kN\n\n"
        f"**Answer:** The maximum reaction at A is {R_max:.4f} kN"
    )

    return question, solution
