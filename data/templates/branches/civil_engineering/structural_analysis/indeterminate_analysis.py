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


# Template 19 (Advanced) — Area A4: Statically Indeterminate Analysis
def template_force_method_continuous_beam():
    """
    Middle-Support Reaction of a Two-Span Continuous Beam (Force Method)

    Scenario:
        A continuous beam rests on supports A, B, C with two equal spans L
        (B at the center). It is statically indeterminate to the first
        degree. Choosing the middle reaction By as the redundant, the
        released structure is a simply supported beam of span 2L, and
        compatibility of deflection at B gives (with delta_B0 the downward
        released-beam deflection and delta_BB the deflection per unit
        UPWARD redundant):

            delta_B0 - By * delta_BB = 0  (deflection at B must vanish)

        The load case changes the released-structure deflection and the
        closed-form result:

            uniform w on both spans:  delta_B0 = 5*w*(2L)^4/(384*EI)
                                      -> By = 5*w*L/4
            point P at each midspan:  delta_B0 = 11*P*L^3/(48*EI)
                                      -> By = 11*P/8

        EI cancels in the compatibility ratio, so no section is needed.

    Difficulty: Advanced
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 10
        (Analysis of Statically Indeterminate Structures by the Force
        Method) — two-span continuous beam with the middle reaction as
        redundant; released-beam deflection formulas per Ch. 8.
    Physical bounds: span L in [3.0, 6.0] m; w in [8, 30] kN/m or P in
        [20, 80] kN per span; By strictly positive and larger than either
        end reaction; end reactions positive.

    Trace integrity (Layer 0, 2026-09-23):
        The end reactions (total load - By)/2 (a 1-dp total minus a 2-dp
        By, halved) are exact at 3 dp and were printed at 2 dp, a half-way
        tie whenever the 2-dp difference has an odd last digit
        ((100.0 - 68.75) / 2 = 15.625, ~44% of draws); they are now bound
        half-up at 3 dp and printed at 3 dp. By = delta_B0/delta_BB is a
        quotient of two rounded displays and exact at no fixed display, so
        a draw whose By lands on a 2-dp half-way tie is resampled rather
        than rounded either way (D-016/D-037).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (load-case branch).
    # Bounded redraw: a draw whose middle reaction lands on a 2-dp half-way
    # tie has no defensible gold answer and is rejected (D-016).
    for _attempt in range(200):
        L = round(random.uniform(3.0, 6.0), 1)
        load_case = random.choice(["uniform", "point"])

        # 2. Core computation — round-then-recompute through the displayed
        # flexibility coefficients (EI carried symbolically; it cancels).
        d_BB = round((2 * L) ** 3 / 48, 3)          # m^3 (times 1/EI), per kN
        if load_case == "uniform":
            w = random.randint(8, 30)
            d_B0 = round(5 * w * (2 * L) ** 4 / 384, 1)  # kN*m^3 (times 1/EI)
            load_text = (f"a uniformly distributed load of {w} kN/m over both "
                         f"spans")
            released_step = (
                f"**Step 2:** Deflection of the released beam at B under the "
                f"real load.\n"
                f"The released structure is a simply supported beam of span "
                f"2L = {2 * L:.1f} m under w = {w} kN/m; its midspan "
                f"deflection is\n"
                f"delta_B0 = 5*w*(2L)^4 / (384*EI) = 5 * {w} * "
                f"({2 * L:.1f})^4 / 384 / EI = {d_B0:.1f}/EI  (kN*m^3 over EI; downward)\n\n"
            )
            total_load = round(2 * w * L, 1)
        else:
            P = random.randint(20, 80)
            d_B0 = round(11 * P * L ** 3 / 48, 1)    # 2*P*(L/2)*(11L^2)/48
            load_text = (f"a concentrated load of {P} kN at the middle of "
                         f"each span")
            released_step = (
                f"**Step 2:** Deflection of the released beam at B under the "
                f"real loads.\n"
                f"The released structure is a simply supported beam of span "
                f"2L = {2 * L:.1f} m carrying {P} kN at x = L/2 = "
                f"{L / 2:.2f} m and at x = 3L/2 = {3 * L / 2:.2f} m. Using "
                f"the standard off-center-load formula "
                f"delta_c = P*b*(3*(2L)^2 - 4*b^2)/(48*EI) with b = L/2 for "
                f"each load and superposing the two equal contributions:\n"
                f"delta_B0 = 2 * {P} * {L / 2:.2f} * (3 * ({2 * L:.1f})^2 - "
                f"4 * ({L / 2:.2f})^2) / 48 / EI = {d_B0:.1f}/EI  "
                f"(kN*m^3 over EI; downward)\n\n"
            )
            total_load = round(2 * P, 1)

        By_exact = d_B0 / d_BB
        if _is_display_tie(By_exact, 2):
            continue                    # no defensible gold answer; redraw
        By = round(By_exact, 2)
        # (total_load - By)/2 is exact at 3 dp; at 2 dp it tied on ~44% of
        # draws (an odd last digit in the 2-dp difference, halved).
        end_R = _hu((total_load - By) / 2, 3)
        break
    else:
        raise RuntimeError("force_method_continuous_beam: no closing sample in 200 draws")

    assert By > 0 and end_R > 0, f"reactions invalid: {By}, {end_R}"
    assert By > end_R, f"middle reaction should dominate: {By} vs {end_R}"

    # 3. Serialize.
    question = (
        f"A continuous beam ABC of constant EI rests on a pin support at "
        f"A, a roller at B, and a roller at C. The two spans are equal: "
        f"AB = BC = L = {L:.1f} m, so B is at the center. The beam "
        f"carries {load_text}. Using the force method with the reaction "
        f"at B as the redundant, determine the vertical reaction at "
        f"support B in kN."
    )

    solution = (
        f"**Given:**\n"
        f"Equal spans: L = {L:.1f} m (total length 2L = {2 * L:.1f} m); "
        f"constant EI\n"
        f"Load: {load_text}\n\n"
        f"**Step 1:** Establish the degree of indeterminacy and release "
        f"the redundant.\n"
        f"Three vertical reactions with two equilibrium equations "
        f"(vertical force and moment) make the beam indeterminate to the "
        f"first degree. Choose By as the redundant and remove support B; "
        f"the released structure is a simply supported beam of span "
        f"2L.\n\n"
        f"{released_step}"
        f"**Step 3:** Flexibility coefficient — deflection at B due to a "
        f"unit upward load at B on the released beam.\n"
        f"delta_BB = (2L)^3 / (48*EI) = ({2 * L:.1f})^3 / 48 / EI "
        f"= {d_BB:.3f}/EI  (m^3 over EI, per kN of redundant)\n\n"
        f"**Step 4:** Apply compatibility at B: the net deflection at "
        f"the support must be zero.\n"
        f"delta_B0 - By * delta_BB = 0\n"
        f"By = delta_B0 / delta_BB = {d_B0:.1f} / {d_BB:.3f} "
        f"= {By:.2f} kN  (EI cancels)\n\n"
        f"**Step 5:** Check equilibrium for the end reactions.\n"
        f"Ay = Cy = (total load - By)/2 = ({total_load:.1f} - {By:.2f}) "
        f"/ 2 = {end_R:.3f} kN, both positive, confirming a consistent "
        f"solution.\n\n"
        f"**Answer:** The vertical reaction at support B is {By:.2f} kN"
    )

    return question, solution


# Template 20 (Advanced) — Area A4: Statically Indeterminate Analysis
def template_slope_deflection_end_moment():
    """
    Fixed-End Moment of a Two-Span Beam by the Slope-Deflection Method

    Scenario:
        A beam ABC has fixed ends at A and C and an interior roller at B,
        with unequal spans L1 and L2 under a uniform load w. The only
        unknown displacement is the joint rotation theta_B. With
        theta_A = theta_C = 0 and no sidesway, the slope-deflection
        equations reduce to a single joint-equilibrium equation in
        X = EI*theta_B:

            M_BA = (4/L1)*X + w*L1^2/12
            M_BC = (4/L2)*X - w*L2^2/12
            M_BA + M_BC = 0  ->  X
            M_AB = (2/L1)*X - w*L1^2/12

        (Clockwise end moments positive; EI cancels in the final moments.)

    Difficulty: Advanced
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 11
        (Displacement Method of Analysis: Slope-Deflection Equations) —
        continuous beam with fixed far ends, single unknown rotation.
    Physical bounds: the shorter span is in [3.0, 5.5] m and the longer
        exceeds it by [0.8, 2.5] m; either may be span AB, so L1 realizes
        values up to 8.0 m (spans must differ — equal spans would make
        theta_B = 0 and the problem degenerate, lesson 17); w sampled in
        [10, 35] kN/m inside a per-sample window that keeps |M_AB| in
        [5, 200] kN*m at both ends.

    Trace integrity (Layer 0, 2026-09-23):
        The fixed-end moments w*L^2/12 were printed at 2 dp. With an
        integer w and a 1-dp span, w*L^2 has 2 dp, so the quotient by 12
        either terminates within 4 dp (when 3 divides 100*w*L^2) or never
        terminates: at 2 dp it sat on a half-way tie on ~7% of draws per
        span (-18 * 5.5^2 / 12 = -45.375), and a 3-dp display would tie on
        the odd terminating cases instead (15 * 3.5^2 / 12 = 15.3125). They
        are therefore bound half-up at 4 dp and printed at 4 dp everywhere,
        as is the joint moment sum rhs, so nothing rounds and no tie can
        arise. X = rhs/coeff and M_AB = (2/L1)*X + FEM_AB are quotients by
        arbitrary rounded divisors, exact at no fixed display, so a draw
        that lands on a half-way tie at its display (3 dp for X, 2 dp for
        M_AB) is resampled rather than rounded either way (D-016/D-037).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize; spans must differ so theta_B != 0.
    # Bounded redraw: a draw whose X or M_AB lands on a half-way tie at its
    # display has no defensible gold answer and is rejected (D-016).
    for _attempt in range(200):
        La = round(random.uniform(3.0, 5.5), 1)
        Lb = round(La + random.uniform(0.8, 2.5), 1)
        if random.random() < 0.5:
            L1, L2 = La, Lb          # shorter span at the fixed end A
        else:
            L1, L2 = Lb, La
        # Per-sample floor on w (lesson 1): the closed form gives
        # |M_AB| = w*|L2*(L2-L1) - 2*L1^2|/24, and the geometry expression is
        # bounded away from zero (min ~4.25 over the sampling space; the
        # cancellation root L2 = 2*L1 is unreachable since L1 >= 3.0 and the
        # spread <= 2.5), so w >= ceil(5.3*24/|expr|) <= 30 keeps |M_AB|
        # above the 5.0 kN*m floor with rounding margin.
        expr = abs(L2 * (L2 - L1) - 2 * L1 ** 2)
        w_lo = max(10, math.ceil(5.3 * 24 / expr))
        # Ceiling too (R1, cycle 1): the long-span branch (expr up to 141.75)
        # with w = 34-35 pushed |M_AB| past the 200 kN*m assert (~1 in 16k
        # draws). w_hi = floor(199*24/expr) >= 33 there, and w_hi/w_lo ~ 37.5
        # so the window is never empty.
        w_hi = min(35, math.floor(199.0 * 24 / expr))
        w = random.randint(w_lo, w_hi)

        # 2. Core computation — round-then-recompute at every step. The
        # fixed-end moments are exact at 4 dp whenever they terminate (a 2-dp
        # w*L^2 over 4*3), so a 4-dp display is tie-free by construction.
        FEM_AB = _hu(-w * L1 ** 2 / 12, 4)
        FEM_BA = _hu(w * L1 ** 2 / 12, 4)
        FEM_BC = _hu(-w * L2 ** 2 / 12, 4)
        FEM_CB = _hu(w * L2 ** 2 / 12, 4)

        coeff = round(4 / L1 + 4 / L2, 4)
        rhs = _hu(-(FEM_BA + FEM_BC), 4)
        X_exact = rhs / coeff
        if _is_display_tie(X_exact, 3):
            continue                    # no defensible gold answer; redraw
        X = round(X_exact, 3)                     # X = EI*theta_B, kN*m^2
        M_exact = 2 / L1 * X + FEM_AB
        if _is_display_tie(M_exact, 2):
            continue
        M_AB = round(M_exact, 2)
        break
    else:
        raise RuntimeError("slope_deflection_end_moment: no closing sample in 200 draws")

    assert abs(FEM_BA + FEM_BC) >= 3.0, "near-degenerate rotation"
    assert 5.0 <= abs(M_AB) <= 200.0, f"end moment out of bounds: {M_AB}"

    # 3. Serialize.
    question = (
        f"A beam ABC of constant EI is fixed at A, rests on a roller at "
        f"B, and is fixed at C. Span AB has length L1 = {L1:.1f} m and "
        f"span BC has length L2 = {L2:.1f} m. A uniformly distributed "
        f"load of w = {w} kN/m acts over both spans. There is no support "
        f"settlement or sidesway. Using the slope-deflection method with "
        f"the standard sign convention (clockwise end moments positive) "
        f"and fixed-end moments FEM = -w*L^2/12 (near end) and "
        f"+w*L^2/12 (far end) for a uniformly loaded span, determine the "
        f"end moment M_AB at the fixed support A in kN*m."
    )

    solution = (
        f"**Given:**\n"
        f"Spans: L1 = {L1:.1f} m (AB), L2 = {L2:.1f} m (BC); w = {w} "
        f"kN/m on both; theta_A = theta_C = 0 (fixed ends), no "
        f"sidesway.\n\n"
        f"**Step 1:** Compute the fixed-end moments for each span.\n"
        f"FEM_AB = -w*L1^2/12 = -{w} * ({L1:.1f})^2 / 12 = "
        f"{FEM_AB:.4f} kN*m\n"
        f"FEM_BA = +w*L1^2/12 = {FEM_BA:.4f} kN*m\n"
        f"FEM_BC = -w*L2^2/12 = -{w} * ({L2:.1f})^2 / 12 = "
        f"{FEM_BC:.4f} kN*m\n"
        f"FEM_CB = +w*L2^2/12 = {FEM_CB:.4f} kN*m\n\n"
        f"**Step 2:** Write the slope-deflection equations with "
        f"theta_A = theta_C = 0, letting X = EI*theta_B.\n"
        f"M_BA = (4/L1)*X + FEM_BA = (4/{L1:.1f})*X + {FEM_BA:.4f}\n"
        f"M_BC = (4/L2)*X + FEM_BC = (4/{L2:.1f})*X + ({FEM_BC:.4f})\n\n"
        f"**Step 3:** Enforce moment equilibrium of joint B.\n"
        f"M_BA + M_BC = 0. Collecting the X terms and the fixed-end "
        f"moments:\n"
        f"(4/{L1:.1f} + 4/{L2:.1f})*X = -({FEM_BA:.4f} + ({FEM_BC:.4f}))\n"
        f"({coeff:.4f})*X = {rhs:.4f}\n"
        f"X = {rhs:.4f} / {coeff:.4f} = {X:.3f} kN*m^2  (EI*theta_B)\n\n"
        f"**Step 4:** Back-substitute into the slope-deflection equation "
        f"for the moment at A.\n"
        f"M_AB = (2/L1)*X + FEM_AB = (2/{L1:.1f}) * {X:.3f} + "
        f"({FEM_AB:.4f}) = {M_AB:.2f} kN*m\n\n"
        f"**Answer:** The end moment at support A is {M_AB:.2f} kN*m"
    )

    return question, solution
