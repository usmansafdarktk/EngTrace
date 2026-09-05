import math
import random


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

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (load-case branch).
    L = round(random.uniform(3.0, 6.0), 1)
    load_case = random.choice(["uniform", "point"])

    # 2. Core computation — round-then-recompute through the displayed
    # flexibility coefficients (EI carried symbolically; it cancels).
    d_BB = round((2 * L) ** 3 / 48, 3)              # m^3 (times 1/EI), per kN
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
        d_B0 = round(11 * P * L ** 3 / 48, 1)        # 2*P*(L/2)*(11L^2)/48
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

    By = round(d_B0 / d_BB, 2)
    end_R = round((total_load - By) / 2, 2)

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
        f"/ 2 = {end_R:.2f} kN, both positive, confirming a consistent "
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

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize; spans must differ so theta_B != 0.
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

    # 2. Core computation — round-then-recompute at every step.
    FEM_AB = round(-w * L1 ** 2 / 12, 2)
    FEM_BA = round(w * L1 ** 2 / 12, 2)
    FEM_BC = round(-w * L2 ** 2 / 12, 2)
    FEM_CB = round(w * L2 ** 2 / 12, 2)

    coeff = round(4 / L1 + 4 / L2, 4)
    rhs = round(-(FEM_BA + FEM_BC), 2)
    X = round(rhs / coeff, 3)                     # X = EI*theta_B, kN*m^2
    M_AB = round(2 / L1 * X + FEM_AB, 2)

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
        f"{FEM_AB:.2f} kN*m\n"
        f"FEM_BA = +w*L1^2/12 = {FEM_BA:.2f} kN*m\n"
        f"FEM_BC = -w*L2^2/12 = -{w} * ({L2:.1f})^2 / 12 = "
        f"{FEM_BC:.2f} kN*m\n"
        f"FEM_CB = +w*L2^2/12 = {FEM_CB:.2f} kN*m\n\n"
        f"**Step 2:** Write the slope-deflection equations with "
        f"theta_A = theta_C = 0, letting X = EI*theta_B.\n"
        f"M_BA = (4/L1)*X + FEM_BA = (4/{L1:.1f})*X + {FEM_BA:.2f}\n"
        f"M_BC = (4/L2)*X + FEM_BC = (4/{L2:.1f})*X + ({FEM_BC:.2f})\n\n"
        f"**Step 3:** Enforce moment equilibrium of joint B.\n"
        f"M_BA + M_BC = 0. Collecting the X terms and the fixed-end "
        f"moments:\n"
        f"(4/{L1:.1f} + 4/{L2:.1f})*X = -({FEM_BA:.2f} + ({FEM_BC:.2f}))\n"
        f"({coeff:.4f})*X = {rhs:.2f}\n"
        f"X = {rhs:.2f} / {coeff:.4f} = {X:.3f} kN*m^2  (EI*theta_B)\n\n"
        f"**Step 4:** Back-substitute into the slope-deflection equation "
        f"for the moment at A.\n"
        f"M_AB = (2/L1)*X + FEM_AB = (2/{L1:.1f}) * {X:.3f} + "
        f"({FEM_AB:.2f}) = {M_AB:.2f} kN*m\n\n"
        f"**Answer:** The end moment at support A is {M_AB:.2f} kN*m"
    )

    return question, solution
