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

    Trace integrity (Layer 0, 2026-09-23):
        The triangular resultant W = w_max*L/2 (an integer times a 1-dp
        span, halved) is exact at 2 dp and was printed at 1 dp, where it
        sat on a half-way tie on about a tenth of triangular draws
        (5 * 9.7 / 2 = 24.25); it is now bound half-up at 2 dp and printed
        at 2 dp everywhere it appears. The uniform resultant w*L is exact
        at 1 dp and is unchanged. The right reaction By = (P*a + W*x_bar)/L
        is a quotient by an arbitrary 1-dp span and exact at no fixed
        display, so a draw whose By lands on a 2-dp half-way tie is
        resampled rather than rounded either way (D-016/D-037).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (load-type branch changes the resultant sub-chain).
    # Bounded redraw: a draw whose right reaction lands on a 2-dp half-way
    # tie has no defensible gold answer and is rejected (D-016).
    for _attempt in range(200):
        L = round(random.uniform(4.0, 10.0), 1)
        P = random.randint(10, 50)
        a = round(random.uniform(0.2 * L, 0.8 * L), 1)
        load_type = random.choice(["uniform", "triangular"])
        if load_type == "uniform":
            w = random.randint(2, 10)
            W = round(w * L, 1)
            W_dp = 1                    # integer w times a 1-dp span: exact
            x_bar = round(L / 2, 2)
            load_text = (f"a uniformly distributed load of {w} kN/m over the "
                         f"entire span")
            step1 = (
                f"**Step 1:** Replace the distributed load by its resultant.\n"
                f"For a uniform load over the full span:\n"
                f"W = w * L = {w} * {L:.1f} = {W:.{W_dp}f} kN, acting at the "
                f"midspan, x_bar = L/2 = {x_bar:.2f} m from A.\n\n"
            )
        else:
            w = random.randint(3, 12)
            # w_max*L/2 is exact at 2 dp; at 1 dp it tied on ~10% of draws.
            W = _hu(w * L / 2, 2)
            W_dp = 2
            x_bar = round(2 * L / 3, 2)
            load_text = (f"a triangularly distributed load that varies from "
                         f"zero at A to {w} kN/m at B")
            step1 = (
                f"**Step 1:** Replace the distributed load by its resultant.\n"
                f"For a triangular load (zero at A, peak at B):\n"
                f"W = w_max * L / 2 = {w} * {L:.1f} / 2 = {W:.{W_dp}f} kN, "
                f"acting at two-thirds of the span from A, x_bar = 2L/3 = "
                f"{x_bar:.2f} m.\n\n"
            )

        # 2. Core computation — round-then-recompute at every step.
        By_exact = (P * a + W * x_bar) / L
        if _is_display_tie(By_exact, 2):
            continue                    # no defensible gold answer; redraw
        By = round(By_exact, 2)
        Ay = round(P + W - By, 2)
        break
    else:
        raise RuntimeError("beam_support_reactions: no closing sample in 200 draws")

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
        f"By = ({P} * {a:.1f} + {W:.{W_dp}f} * {x_bar:.2f}) / {L:.1f} "
        f"= {By:.2f} kN\n\n"
        f"**Step 3:** Apply vertical force equilibrium to find the "
        f"reaction at A.\n"
        f"Sum(F_y) = 0: Ay = P + W - By = {P} + {W:.{W_dp}f} - {By:.2f} "
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

    Trace integrity (Layer 0, 2026-09-23):
        F_AB = Ay / sin(theta) is a quotient by a 4-dp sine, exact at no
        fixed display; on the 3-4-5 geometries sin(theta) = 0.8000 exactly
        and F_AB = 1.25 * Ay lands on a 2-dp half-way tie whenever P is odd
        (32.5 / 0.8000 = 40.625). Such a draw, and the draw whose sine
        h / L_AB itself lands on a 4-dp tie, is resampled rather than
        rounded either way (D-016/D-037); the answer's 2-dp display is
        unchanged.

    Screen pass 1 (2026-09-23):
        A judge read the gold as resting on an unstated L-to-3-dp /
        sin-to-4-dp rounding path that can differ from the exact force.
        The path is the stated round-then-recompute convention: every
        consumed operand is the value printed, so each line closes from
        the line above. Measured at 500 seeds, the 2-dp gold differs from
        the 2-dp rounding of the exact force on 22.8% of draws, by one
        unit in the last place (two on a single draw; worst relative gap
        4.4e-4). Step 2 now states the path in one sentence; no number,
        draw or step changes. The same judge's "3.0 / 4.243 = 0.7070" is
        not a defect: the exact quotient is 0.7070469..., which rounds to
        0.7070, and the line closes.

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize with the inclination kept in a sensible window.
    # Bounded redraw: a draw whose diagonal force (or its sine) lands on a
    # half-way tie at its display has no defensible gold answer and is
    # rejected (D-016).
    for _attempt in range(200):
        b = round(random.uniform(2.0, 4.0), 1)
        h = round(random.uniform(0.47 * b, 1.7 * b), 1)   # theta ~ 25-60 deg
        P = random.randint(20, 80)

        # 2. Core computation — round-then-recompute at every step.
        L_ab = round(math.sqrt(h ** 2 + b ** 2), 3)
        sin_exact = h / L_ab
        sin_t = round(sin_exact, 4)
        Ay = round(P / 2, 1)
        F_exact = Ay / sin_t
        if _is_display_tie(sin_exact, 4) or _is_display_tie(F_exact, 2):
            continue                    # no defensible gold answer; redraw
        F_ab = round(F_exact, 2)
        break
    else:
        raise RuntimeError("truss_method_of_joints: no closing sample in 200 draws")

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
        f"Lengths are rounded to 3 dp and sines to 4 dp; each later step "
        f"uses the rounded value as printed.\n"
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

    Trace integrity (Layer 0, 2026-09-23):
        The moment M = Ay*x - w*x^2/2 [- P*(x - a)], with a 2-dp reaction
        and a 1-dp section, is exact at 3 dp and was printed at 2 dp, where
        it sat on a half-way tie on ~9% of draws (26.25 * 2.9 - 4 *
        2.9^2 / 2 = 59.305). M is the answer, so on the owner's decision of
        2026-09-23 (D-044) its display was lengthened to 3 dp rather than
        the draw resampled: M is bound half-up at 3 dp and printed at 3 dp
        in Step 3 and in the answer, which removes the tie by
        construction. By = (P*a + w*L^2/2)/L is a quotient by an
        arbitrary 1-dp span, exact at no fixed display, so a draw whose By
        lands on a 2-dp tie is resampled rather than rounded either way
        (D-016/D-037).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (section-side branch).
    # Bounded redraw: a draw whose reaction By lands on a 2-dp half-way tie
    # has no defensible gold answer and is rejected (D-016).
    for _attempt in range(200):
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
        By_exact = (P * a + w * L * L / 2) / L
        By = round(By_exact, 2)
        Ay = round(P + w * L - By, 2)
        if x > a:
            M_exact = Ay * x - w * x ** 2 / 2 - P * (x - a)
        else:
            M_exact = Ay * x - w * x ** 2 / 2
        if _is_display_tie(By_exact, 2):
            continue                    # no defensible gold answer; redraw
        # M is exact at 3 dp (2-dp reaction times 1-dp section, integer
        # loads); bound half-up and printed at 3 dp (D-044).
        M = _hu(M_exact, 3)
        break
    else:
        raise RuntimeError("beam_internal_moment: no closing sample in 200 draws")

    if x > a:
        cut_terms = (
            f"the reaction Ay, the distributed load over the length x, "
            f"and the point load P (since x = {x:.1f} m > a = {a:.1f} m, "
            f"the point load acts on the retained segment)")
        moment_line = (
            f"M = Ay * x - w * x^2 / 2 - P * (x - a)\n"
            f"M = {Ay:.2f} * {x:.1f} - {w} * ({x:.1f})^2 / 2 - {P} * "
            f"({x:.1f} - {a:.1f}) = {M:.3f} kN*m")
    else:
        cut_terms = (
            f"the reaction Ay and the distributed load over the length x "
            f"only (the point load acts beyond the section, since "
            f"x = {x:.1f} m < a = {a:.1f} m, so it does not load the "
            f"retained segment)")
        moment_line = (
            f"M = Ay * x - w * x^2 / 2\n"
            f"M = {Ay:.2f} * {x:.1f} - {w} * ({x:.1f})^2 / 2 "
            f"= {M:.3f} kN*m")

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
        f"{M:.3f} kN*m"
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

    Trace integrity (Layer 0, 2026-09-23):
        Every member force is a quotient by the arbitrary 1-dp height h
        (Ay*d/h, Ay*2d/h, Ay*L_BC/h), exact at no fixed display, and lands
        on a 2-dp half-way tie on ~3% of draws, concentrated on the heights
        whose 10*h has only the factors 2 and 5 (33.5 * 2.0 / 1.6 =
        41.875). The force is the answer, so its 2-dp display is kept and
        such a draw is resampled rather than rounded either way
        (D-016/D-037).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (member branch selects the equilibrium equation).
    # Bounded redraw: a draw whose member force lands on a 2-dp half-way
    # tie has no defensible gold answer and is rejected (D-016).
    for _attempt in range(200):
        d = round(random.uniform(2.0, 4.0), 1)
        h = round(random.uniform(0.6 * d, 1.2 * d), 1)
        P = random.randint(20, 80)
        member = random.choice(["AC", "BD", "BC"])

        # 2. Core computation — round-then-recompute at every step.
        Ay = round(P / 2, 1)
        if member == "AC":
            F_exact = Ay * d / h
        elif member == "BD":
            F_exact = Ay * 2 * d / h
        else:
            L_bc = round(math.sqrt(d ** 2 + h ** 2), 3)
            F_exact = Ay * L_bc / h
        if _is_display_tie(F_exact, 2):
            continue                    # no defensible gold answer; redraw
        F = round(F_exact, 2)
        break
    else:
        raise RuntimeError("truss_method_of_sections: no closing sample in 200 draws")

    if member == "AC":
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
