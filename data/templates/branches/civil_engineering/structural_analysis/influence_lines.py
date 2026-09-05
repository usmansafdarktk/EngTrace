import random


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

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize (heavier-axle identity is the positioning branch).
    L = round(random.uniform(8.0, 16.0), 1)
    P1 = random.randint(40, 160)                  # front axle
    P2 = random.choice([p for p in range(40, 161)
                        if abs(p - P1) >= 10])    # rear axle
    s = round(random.uniform(2.0, min(4.5, L / 3 - 0.1)), 1)

    P_heavy, P_light = (P1, P2) if P1 > P2 else (P2, P1)
    heavy_name = "front" if P1 > P2 else "rear"

    # 2. Core computation — round-then-recompute at every step.
    y_light = round((L - s) / L, 4)
    R_max = round(P_heavy * 1.0 + P_light * y_light, 2)

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
        f"= {R_max:.2f} kN\n\n"
        f"**Answer:** The maximum reaction at A is {R_max:.2f} kN"
    )

    return question, solution
