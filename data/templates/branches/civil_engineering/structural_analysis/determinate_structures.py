import math
import random


# Template 11 (Easy) — Area A1: Analysis of Determinate Structures
def template_beam_support_reactions():
    """
    Support Reactions of a Simply Supported Beam

    Scenario:
        A simply supported beam carries a point load plus a distributed
        load over the full span. The distributed load is either uniform or
        triangular (zero at the left support, peak at the right), which
        changes both the resultant magnitude and its line of action:

            uniform:    W = w*L   at x = L/2
            triangular: W = w_max*L/2   at x = 2L/3

        Moment equilibrium about A gives the right reaction; vertical
        force equilibrium gives the left reaction (the requested answer).

    Difficulty: Easy
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Section 2.4
        (Equations of Equilibrium); distributed-load resultants per Ch. 2
        idealized-load treatment.
    Physical bounds: span L in [4.0, 10.0] m; point load P in [10, 50] kN
        at a in [0.2L, 0.8L]; distributed intensity in [2, 10] kN/m
        (uniform) or peak [3, 12] kN/m (triangular); both reactions
        strictly positive and their sum equals the total load within
        rounding.

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (load-type branch changes the resultant sub-chain).
    L = round(random.uniform(4.0, 10.0), 1)
    P = random.randint(10, 50)
    a = round(random.uniform(0.2 * L, 0.8 * L), 1)
    load_type = random.choice(["uniform", "triangular"])
    if load_type == "uniform":
        w = random.randint(2, 10)
        W = round(w * L, 1)
        x_bar = round(L / 2, 2)
        load_text = (f"a uniformly distributed load of {w} kN/m over the "
                     f"entire span")
        step1 = (
            f"**Step 1:** Replace the distributed load by its resultant.\n"
            f"For a uniform load over the full span:\n"
            f"W = w * L = {w} * {L:.1f} = {W:.1f} kN, acting at the "
            f"midspan, x_bar = L/2 = {x_bar:.2f} m from A.\n\n"
        )
    else:
        w = random.randint(3, 12)
        W = round(0.5 * w * L, 1)
        x_bar = round(2 * L / 3, 2)
        load_text = (f"a triangularly distributed load that varies from "
                     f"zero at A to {w} kN/m at B")
        step1 = (
            f"**Step 1:** Replace the distributed load by its resultant.\n"
            f"For a triangular load (zero at A, peak at B):\n"
            f"W = w_max * L / 2 = {w} * {L:.1f} / 2 = {W:.1f} kN, acting "
            f"at two-thirds of the span from A, x_bar = 2L/3 = "
            f"{x_bar:.2f} m.\n\n"
        )

    # 2. Core computation — round-then-recompute at every step.
    By = round((P * a + W * x_bar) / L, 2)
    Ay = round(P + W - By, 2)

    assert Ay > 0 and By > 0, f"reaction not positive: {Ay}, {By}"
    assert abs((Ay + By) - (P + W)) < 0.02, "equilibrium check failed"

    # 3. Serialize.
    question = (
        f"A simply supported beam AB has a span of L = {L:.1f} m, with a "
        f"pin support at A and a roller support at B. It carries a "
        f"vertical downward concentrated load of {P} kN at a distance of "
        f"{a:.1f} m from A, together with {load_text} (also acting "
        f"downward). Determine the vertical reaction at support A in kN."
    )

    solution = (
        f"**Given:**\n"
        f"Span (L): {L:.1f} m\n"
        f"Point load (P): {P} kN at a = {a:.1f} m from A\n"
        f"Distributed load: {load_text}\n\n"
        f"{step1}"
        f"**Step 2:** Take moments about A to find the reaction at B.\n"
        f"Sum(M_A) = 0: By * L = P * a + W * x_bar\n"
        f"By = ({P} * {a:.1f} + {W:.1f} * {x_bar:.2f}) / {L:.1f} "
        f"= {By:.2f} kN\n\n"
        f"**Step 3:** Apply vertical force equilibrium to find the "
        f"reaction at A.\n"
        f"Sum(F_y) = 0: Ay = P + W - By = {P} + {W:.1f} - {By:.2f} "
        f"= {Ay:.2f} kN\n\n"
        f"**Answer:** The vertical reaction at support A is {Ay:.2f} kN"
    )

    return question, solution


# Template 12 (Easy) — Area A1: Analysis of Determinate Structures
def template_truss_method_of_joints():
    """
    Diagonal Member Force by the Method of Joints

    Scenario:
        A symmetric triangular truss carries a single load at its apex.
        Symmetry gives the support reactions; equilibrium of the support
        joint then yields the diagonal force from the member's inclination:

            sin(theta) = h / sqrt(h^2 + b^2)
            F_diag = Ay / sin(theta)   (compression)

    Difficulty: Easy
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Section 3.3
        (The Method of Joints).
    Physical bounds: half-span b in [2.0, 4.0] m; height h chosen so the
        diagonal's inclination stays in [24, 61] degrees (1-dp rounding of
        h can nudge the nominal 25-60 window by ~1 degree); apex load P in
        [20, 80] kN; the diagonal is always in compression.

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize with the inclination kept in a sensible window.
    b = round(random.uniform(2.0, 4.0), 1)
    h = round(random.uniform(0.47 * b, 1.7 * b), 1)   # theta ~ 25-60 deg
    P = random.randint(20, 80)

    # 2. Core computation — round-then-recompute at every step.
    L_ab = round(math.sqrt(h ** 2 + b ** 2), 3)
    sin_t = round(h / L_ab, 4)
    Ay = round(P / 2, 1)
    F_ab = round(Ay / sin_t, 2)

    theta_deg = math.degrees(math.asin(min(1.0, sin_t)))
    assert 24.0 <= theta_deg <= 61.0, f"inclination out of window: {theta_deg}"
    assert F_ab > Ay, "diagonal force must exceed the reaction"

    # 3. Serialize.
    question = (
        f"A symmetric triangular truss has supports A (pin) and C "
        f"(roller) at the same level, a horizontal distance 2b = "
        f"{2 * b:.1f} m apart, and an apex joint B at height h = "
        f"{h:.1f} m above the midpoint of AC. The three members are the "
        f"diagonals AB and BC and the bottom chord AC. A single vertical "
        f"downward load of P = {P} kN acts at the apex B. Using the "
        f"method of joints, determine the magnitude of the force in "
        f"diagonal member AB in kN."
    )

    solution = (
        f"**Given:**\n"
        f"Half-span (b): {b:.1f} m, apex height (h): {h:.1f} m\n"
        f"Apex load (P): {P} kN\n\n"
        f"**Step 1:** Find the support reactions.\n"
        f"By symmetry of the truss and load, each support carries half "
        f"the load:\n"
        f"Ay = P / 2 = {P} / 2 = {Ay:.1f} kN\n\n"
        f"**Step 2:** Establish the geometry of diagonal AB.\n"
        f"Let theta be the inclination of AB to the horizontal chord AC.\n"
        f"Member length: L_AB = sqrt(h^2 + b^2) = sqrt({h:.1f}^2 + "
        f"{b:.1f}^2) = {L_ab:.3f} m\n"
        f"sin(theta) = h / L_AB = {h:.1f} / {L_ab:.3f} = {sin_t:.4f}\n\n"
        f"**Step 3:** Apply vertical equilibrium at joint A.\n"
        f"Only AB has a vertical component at A (AC is horizontal). To "
        f"balance the upward reaction Ay, AB's force on the joint must "
        f"act downward along the member, which means the member pushes "
        f"on the joint — AB is in compression:\n"
        f"F_AB = Ay / sin(theta) = {Ay:.1f} / {sin_t:.4f} = {F_ab:.2f} kN "
        f"(compression)\n\n"
        f"**Answer:** The force in member AB is {F_ab:.2f} kN"
    )

    return question, solution


# Template 13 (Easy) — Area A1: Analysis of Determinate Structures
def template_beam_internal_moment():
    """
    Internal Bending Moment at a Beam Section

    Scenario:
        A simply supported beam carries a full-span uniform load and a
        point load. The requested section lies either left or right of the
        point load — a branch that changes which forces cross the cut and
        therefore the moment expression:

            x < a:  M = Ay*x - w*x^2/2
            x > a:  M = Ay*x - w*x^2/2 - P*(x - a)

    Difficulty: Easy
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 4
        (Internal Loadings Developed in Structural Members) with reactions
        per Section 2.4.
    Physical bounds: span L in [6.0, 12.0] m; P in [15, 60] kN at
        a in [0.3L, 0.7L]; w in [2, 8] kN/m; section location kept at
        least 0.5 m from both supports and from the load point.

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (section-side branch).
    L = round(random.uniform(6.0, 12.0), 1)
    P = random.randint(15, 60)
    a = round(random.uniform(0.3 * L, 0.7 * L), 1)
    w = random.randint(2, 8)
    side = random.choice(["left", "right"])
    if side == "left":
        x = round(random.uniform(0.15 * L, a - 0.5), 1)
    else:
        x = round(random.uniform(a + 0.5, 0.85 * L), 1)

    # 2. Core computation — round-then-recompute at every step.
    By = round((P * a + w * L * L / 2) / L, 2)
    Ay = round(P + w * L - By, 2)
    if x > a:
        M = round(Ay * x - w * x ** 2 / 2 - P * (x - a), 2)
        cut_terms = (
            f"the reaction Ay, the distributed load over the length x, "
            f"and the point load P (since x = {x:.1f} m > a = {a:.1f} m, "
            f"the point load acts on the retained segment)")
        moment_line = (
            f"M = Ay * x - w * x^2 / 2 - P * (x - a)\n"
            f"M = {Ay:.2f} * {x:.1f} - {w} * ({x:.1f})^2 / 2 - {P} * "
            f"({x:.1f} - {a:.1f}) = {M:.2f} kN*m")
    else:
        M = round(Ay * x - w * x ** 2 / 2, 2)
        cut_terms = (
            f"the reaction Ay and the distributed load over the length x "
            f"only (the point load acts beyond the section, since "
            f"x = {x:.1f} m < a = {a:.1f} m, so it does not load the "
            f"retained segment)")
        moment_line = (
            f"M = Ay * x - w * x^2 / 2\n"
            f"M = {Ay:.2f} * {x:.1f} - {w} * ({x:.1f})^2 / 2 "
            f"= {M:.2f} kN*m")

    assert 0.5 <= x <= L - 0.5 and abs(x - a) >= 0.45, (
        f"section location invalid: {x}")
    assert M > 0, f"sagging moment expected on a simple beam: {M}"

    # 3. Serialize.
    question = (
        f"A simply supported beam AB (pin at A, roller at B) has a span "
        f"of L = {L:.1f} m. It carries a uniformly distributed load of "
        f"{w} kN/m over the entire span and a concentrated load of "
        f"{P} kN at a = {a:.1f} m from A; all loads act vertically "
        f"downward. Determine the internal bending moment (sagging "
        f"positive), in kN*m, at the section located x = {x:.1f} m from "
        f"support A."
    )

    solution = (
        f"**Given:**\n"
        f"Span (L): {L:.1f} m; UDL (w): {w} kN/m over the full span\n"
        f"Point load (P): {P} kN at a = {a:.1f} m from A\n"
        f"Section location: x = {x:.1f} m from A\n\n"
        f"**Step 1:** Find the support reactions.\n"
        f"Sum(M_A) = 0: By = (P * a + w * L * L/2) / L = ({P} * {a:.1f} "
        f"+ {w} * {L:.1f} * {L / 2:.2f}) / {L:.1f} = {By:.2f} kN\n"
        f"Sum(F_y) = 0: Ay = P + w * L - By = {P} + {w} * {L:.1f} - "
        f"{By:.2f} = {Ay:.2f} kN\n\n"
        f"**Step 2:** Cut the beam at x = {x:.1f} m and take the left "
        f"free body.\n"
        f"The external forces acting on the retained left segment are "
        f"{cut_terms}.\n\n"
        f"**Step 3:** Sum moments about the cut to find the internal "
        f"moment (sagging positive).\n"
        f"{moment_line}\n\n"
        f"**Answer:** The internal bending moment at the section is "
        f"{M:.2f} kN*m"
    )

    return question, solution


# Template 14 (Intermediate) — Area A1: Analysis of Determinate Structures
def template_truss_method_of_sections():
    """
    Member Force in a Warren Truss by the Method of Sections

    Scenario:
        A five-joint Warren truss (bottom joints A, C, E; top joints B, D)
        carries a single load at the central bottom joint. A vertical cut
        between B and C crosses the top chord BD, the diagonal BC, and the
        bottom chord AC. The requested member selects the equilibrium
        equation:

            AC: moments about B  ->  F_AC = Ay * d / h        (tension)
            BD: moments about C  ->  F_BD = Ay * 2d / h       (compression)
            BC: vertical forces  ->  F_BC = Ay * L_BC / h     (tension)

    Difficulty: Intermediate
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Section 3.5
        (The Method of Sections).
    Physical bounds: panel length d in [2.0, 4.0] m; height h in
        [0.6d, 1.2d] (diagonal inclination ~31-50 deg); load P in
        [20, 80] kN; requested member sampled from {AC, BD, BC}.

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (member branch selects the equilibrium equation).
    d = round(random.uniform(2.0, 4.0), 1)
    h = round(random.uniform(0.6 * d, 1.2 * d), 1)
    P = random.randint(20, 80)
    member = random.choice(["AC", "BD", "BC"])

    # 2. Core computation — round-then-recompute at every step.
    Ay = round(P / 2, 1)
    if member == "AC":
        F = round(Ay * d / h, 2)
        sense = "tension"
        method_step = (
            f"**Step 3:** Take moments about joint B (at x = d, height h) "
            f"for the left free body.\n"
            f"Only Ay and F_AC have moment about B (BD and BC pass "
            f"through B):\n"
            f"F_AC * h = Ay * d\n"
            f"F_AC = Ay * d / h = {Ay:.1f} * {d:.1f} / {h:.1f} = "
            f"{F:.2f} kN (tension)")
    elif member == "BD":
        F = round(Ay * 2 * d / h, 2)
        sense = "compression"
        method_step = (
            f"**Step 3:** Take moments about joint C (at x = 2d on the "
            f"base) for the left free body.\n"
            f"Only Ay and F_BD have moment about C (BC and AC pass "
            f"through C):\n"
            f"F_BD * h = Ay * 2d\n"
            f"F_BD = Ay * 2d / h = {Ay:.1f} * {2 * d:.1f} / {h:.1f} = "
            f"{F:.2f} kN (compression)")
    else:
        L_bc = round(math.sqrt(d ** 2 + h ** 2), 3)
        F = round(Ay * L_bc / h, 2)
        sense = "tension"
        method_step = (
            f"**Step 3:** Sum vertical forces on the left free body.\n"
            f"Chords AC and BD are horizontal, so the diagonal BC must "
            f"carry the panel shear. Its length is L_BC = sqrt(d^2 + "
            f"h^2) = sqrt({d:.1f}^2 + {h:.1f}^2) = {L_bc:.3f} m, so its "
            f"vertical component fraction is h / L_BC.\n"
            f"F_BC * (h / L_BC) = Ay\n"
            f"F_BC = Ay * L_BC / h = {Ay:.1f} * {L_bc:.3f} / {h:.1f} = "
            f"{F:.2f} kN (tension)")

    theta = math.degrees(math.atan(h / d))
    assert 30.0 <= theta <= 51.0, f"diagonal inclination out of window: {theta}"
    assert F > 0, f"member force must be positive: {F}"

    # 3. Serialize.
    question = (
        f"A Warren truss has bottom-chord joints A, C, and E at "
        f"x = 0, {2 * d:.1f} m, and {4 * d:.1f} m along the base, and "
        f"top-chord joints B and D at x = {d:.1f} m and {3 * d:.1f} m, "
        f"both at height h = {h:.1f} m. The members are AB, BC, CD, DE, "
        f"the bottom chords AC and CE, and the top chord BD. The truss "
        f"has a pin support at A, a roller support at E, and carries a "
        f"single vertical downward load of P = {P} kN at joint C. Using "
        f"the method of sections with a vertical cut between joints B "
        f"and C (crossing members BD, BC, and AC), determine the "
        f"magnitude of the force in member {member} in kN."
    )

    solution = (
        f"**Given:**\n"
        f"Panel geometry: bottom joints at 0, {2 * d:.1f}, {4 * d:.1f} m; "
        f"top joints at {d:.1f}, {3 * d:.1f} m; height h = {h:.1f} m\n"
        f"Horizontal joint spacing: d = {d:.1f} m (so the bottom joints "
        f"sit at 0, 2d, 4d and the top joints at d, 3d)\n"
        f"Load (P): {P} kN downward at joint C\n"
        f"Requested member: {member}\n\n"
        f"**Step 1:** Find the support reactions.\n"
        f"The load at C is at midspan, so by symmetry:\n"
        f"Ay = Ey = P / 2 = {P} / 2 = {Ay:.1f} kN\n\n"
        f"**Step 2:** Cut the truss between B and C.\n"
        f"The cut crosses the top chord BD, the diagonal BC, and the "
        f"bottom chord AC. Consider the left free body, which carries "
        f"the reaction Ay at A and the three unknown member forces.\n\n"
        f"{method_step}\n\n"
        f"**Answer:** The force in member {member} is {F:.2f} kN "
        f"({sense})"
    )

    return question, solution
