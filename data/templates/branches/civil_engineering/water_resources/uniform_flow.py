import math
import random

from pilot.templates.branches.civil_engineering.constants import (
    GRAVITY_M_S2,
    MANNINGS_N_CHANNELS,
)

# Rigid linings only (R1, cycle 1): the unlined-earth entry drew erosive
# velocities (V up to ~3.9 m/s vs ~1.8 m/s permissible for earth) and
# implied vertical earth walls in rectangular sections, contradicting the
# steady-uniform-flow premise; the conflated "shot concrete / earth
# channel" label is likewise excluded. Values remain the HDS-4 Table B.2
# entries in constants.
_RIGID_LININGS = {
    "very smooth concrete": MANNINGS_N_CHANNELS["very smooth concrete"],
    "smooth concrete": MANNINGS_N_CHANNELS["smooth concrete"],
    "ordinary concrete lining": MANNINGS_N_CHANNELS["ordinary concrete lining"],
    "wood": MANNINGS_N_CHANNELS["wood"],
    "vitrified clay": MANNINGS_N_CHANNELS["vitrified clay"],
}


def _lining_phrase(lining):
    """Natural-language lining phrase: strips the HDS-4 key's own 'lining'
    suffix (R3, cycle 2: 'a ordinary concrete lining lining') and picks
    the right article."""
    display = lining.replace(" lining", "")
    article = "an" if display[0] in "aeiou" else "a"
    return f"{article} {display} lining"


def _froude_capped_slope(n, A, P, T, q_cap_flow, s_floor=0.0005, s_ceil=0.004):
    """Per-sample slope window (lessons 1/14): cap S so the uniform flow
    stays clearly subcritical (Fr <= ~0.88) AND below the discharge cap.
    Returns a rounded slope sample."""
    g = GRAVITY_M_S2
    R = A / P
    D = A / T
    v_lim = 0.88 * math.sqrt(g * D)
    s_fr = 0.94 * (v_lim * n) ** 2 / R ** (4.0 / 3.0)
    s_q = 0.96 * (q_cap_flow * n / (A * R ** (2.0 / 3.0))) ** 2
    s_max = min(s_ceil, s_fr, s_q)
    assert s_max > s_floor, f"empty slope window: {s_max}"
    return round(random.uniform(s_floor, s_max), 4)


# Template 21 (Easy) — Area C1: Uniform Open-Channel Flow
def template_manning_rectangular_discharge():
    """
    Discharge of a Rectangular Channel by Manning's Equation

    Scenario:
        A rigid-lined rectangular channel of known width, flow depth,
        slope, and lining carries uniform flow. The chain is the standard
        Manning computation:

            A = b*y;  P = b + 2y;  R = A/P
            Q = (1/n) * A * R^(2/3) * S^(1/2)

    Difficulty: Easy
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 4 (uniform
        flow; Manning's equation, SI form); n values per
        constants.MANNINGS_N_CHANNELS (FHWA HDS-4 Table B.2, on-disk),
        rigid linings only.
    Physical bounds: width b in [2.0, 6.0] m; depth y in [0.8, 2.5] m;
        slope sampled inside a per-sample window that keeps the uniform
        flow subcritical (Fr <= ~0.88) and the discharge <= ~58 m^3/s;
        discharge asserted in [0.9, 60] m^3/s.

    Returns:
        tuple: (question, solution)
    """
    lining = random.choice(list(_RIGID_LININGS.keys()))
    n = _RIGID_LININGS[lining]
    b = round(random.uniform(2.0, 6.0), 1)
    y = round(random.uniform(0.8, 2.5), 1)
    S = _froude_capped_slope(n, b * y, b + 2 * y, b, 58.0)

    A = round(b * y, 3)
    P = round(b + 2 * y, 2)
    R = round(A / P, 4)
    R23 = round(R ** (2.0 / 3.0), 4)
    sqS = round(math.sqrt(S), 5)
    Q = round((1.0 / n) * A * R23 * sqS, 2)

    Fr = (Q / A) / math.sqrt(GRAVITY_M_S2 * y)
    assert 0.9 <= Q <= 60.0, f"discharge out of bounds: {Q}"
    assert Fr <= 0.92, f"uniform flow not subcritical: Fr = {Fr}"

    question = (
        f"A rectangular open channel is {b:.1f} m wide and has "
        f"{_lining_phrase(lining)} (Manning's n = {n}). The channel carries "
        f"uniform flow at a depth of {y:.1f} m on a longitudinal slope "
        f"of S = {S}. Using Manning's equation in SI units, determine "
        f"the discharge in m^3/s."
    )

    solution = (
        f"**Given:**\n"
        f"Width (b): {b:.1f} m; depth (y): {y:.1f} m\n"
        f"Slope (S): {S}; lining: {lining}, n = {n}\n\n"
        f"**Step 1:** Compute the flow area and wetted perimeter.\n"
        f"A = b * y = {b:.1f} * {y:.1f} = {A:.3f} m^2\n"
        f"P = b + 2y = {b:.1f} + 2 * {y:.1f} = {P:.2f} m\n\n"
        f"**Step 2:** Compute the hydraulic radius.\n"
        f"R = A / P = {A:.3f} / {P:.2f} = {R:.4f} m\n\n"
        f"**Step 3:** Evaluate the Manning terms.\n"
        f"R^(2/3) = ({R:.4f})^(2/3) = {R23:.4f}\n"
        f"S^(1/2) = ({S})^(1/2) = {sqS:.5f}\n\n"
        f"**Step 4:** Apply Manning's equation.\n"
        f"Q = (1/n) * A * R^(2/3) * S^(1/2) "
        f"= (1/{n}) * {A:.3f} * {R23:.4f} * {sqS:.5f} = {Q:.2f} m^3/s\n\n"
        f"**Answer:** The discharge is {Q:.2f} m^3/s"
    )

    return question, solution


# Template 22 (Easy) — Area C1: Uniform Open-Channel Flow
def template_manning_trapezoidal_velocity():
    """
    Mean Velocity in a Trapezoidal Channel by Manning's Equation

    Scenario:
        A rigid-lined trapezoidal channel with side slopes z:1 carries
        uniform flow. The velocity form of Manning's equation terminates
        the chain (distinct from the discharge computation of the
        rectangular template):

            A = (b + z*y)*y;  P = b + 2*y*sqrt(1 + z^2);  R = A/P
            V = (1/n) * R^(2/3) * S^(1/2);  check Q = V*A

    Difficulty: Easy
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 4 (uniform
        flow; velocity form of Manning's equation; trapezoidal geometry);
        n per constants (HDS-4, on-disk), rigid linings only.
    Physical bounds: bottom width b in [1.5, 5.0] m; depth y in
        [0.8, 2.2] m; side slope z in {1.5, 2.0, 2.5, 3.0}; slope sampled
        inside a per-sample subcritical window; velocity asserted in
        [0.4, 3.5] m/s.

    Returns:
        tuple: (question, solution)
    """
    lining = random.choice(list(_RIGID_LININGS.keys()))
    n = _RIGID_LININGS[lining]
    b = round(random.uniform(1.5, 5.0), 1)
    y = round(random.uniform(0.8, 2.2), 1)
    z = random.choice([1.5, 2.0, 2.5, 3.0])
    A_t = (b + z * y) * y
    P_t = b + 2 * y * math.sqrt(1 + z ** 2)
    T_t = b + 2 * z * y
    S = _froude_capped_slope(n, A_t, P_t, T_t, 85.0)

    A = round(A_t, 3)
    root = round(math.sqrt(1 + z ** 2), 4)
    P = round(b + 2 * y * root, 3)
    R = round(A / P, 4)
    R23 = round(R ** (2.0 / 3.0), 4)
    sqS = round(math.sqrt(S), 5)
    V = round((1.0 / n) * R23 * sqS, 3)
    Q_check = round(V * A, 2)

    assert 0.4 <= V <= 3.5, f"velocity out of bounds: {V}"

    question = (
        f"A trapezoidal open channel has a bottom width of {b:.1f} m and "
        f"side slopes of {z}H:1V, with {_lining_phrase(lining)} (Manning's "
        f"n = {n}). It carries uniform flow at a depth of {y:.1f} m on a "
        f"slope of S = {S}. Using the velocity form of Manning's "
        f"equation in SI units, determine the mean flow velocity in m/s."
    )

    solution = (
        f"**Given:**\n"
        f"Bottom width (b): {b:.1f} m; depth (y): {y:.1f} m; side slope "
        f"z = {z}\n"
        f"Slope (S): {S}; lining: {lining}, n = {n}\n\n"
        f"**Step 1:** Compute the flow area.\n"
        f"A = (b + z*y) * y = ({b:.1f} + {z} * {y:.1f}) * {y:.1f} "
        f"= {A:.3f} m^2\n\n"
        f"**Step 2:** Compute the wetted perimeter and hydraulic "
        f"radius.\n"
        f"sqrt(1 + z^2) = sqrt(1 + {z}^2) = {root:.4f}\n"
        f"P = b + 2*y*sqrt(1 + z^2) = {b:.1f} + 2 * {y:.1f} * "
        f"{root:.4f} = {P:.3f} m\n"
        f"R = A / P = {A:.3f} / {P:.3f} = {R:.4f} m\n\n"
        f"**Step 3:** Apply the velocity form of Manning's equation.\n"
        f"R^(2/3) = ({R:.4f})^(2/3) = {R23:.4f}; "
        f"S^(1/2) = ({S})^(1/2) = {sqS:.5f}\n"
        f"V = (1/n) * R^(2/3) * S^(1/2) "
        f"= (1/{n}) * {R23:.4f} * {sqS:.5f} = {V:.3f} m/s\n\n"
        f"**Step 4:** Check the corresponding discharge.\n"
        f"Q = V * A = {V:.3f} * {A:.3f} = {Q_check:.2f} m^3/s, a "
        f"discharge this channel carries at the given depth.\n\n"
        f"**Answer:** The mean flow velocity is {V:.3f} m/s"
    )

    return question, solution


# Template 23 (Intermediate) — Area C1: Uniform Open-Channel Flow
def template_best_hydraulic_rectangular_section():
    """
    Depth of the Best Hydraulic Rectangular Section

    Scenario:
        For a given discharge, slope, and lining, the most efficient
        (best hydraulic) rectangular section has b = 2y and R = y/2.
        Substituting into Manning's equation gives a closed form for the
        required depth:

            Q = (1/n) * 2y^2 * (y/2)^(2/3) * S^(1/2)
              = 2^(1/3) * y^(8/3) * S^(1/2) / n
            y = [Q*n / (2^(1/3) * S^(1/2))]^(3/8)

    Difficulty: Intermediate
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 4 (uniform
        flow; most efficient/best hydraulic sections); n per constants
        (HDS-4, on-disk), rigid linings with n >= 0.012.
    Physical bounds: design discharge Q in [5, 40] m^3/s; S in
        [0.0005, 0.0018] (keeps the resulting flow subcritical, Fr
        asserted <= 0.90); resulting depth asserted in [0.9, 3.5] m; the
        trace closes with a Manning back-check within ~1.5%.

    Returns:
        tuple: (question, solution)
    """
    rigid_no_glassy = {k: v for k, v in _RIGID_LININGS.items() if v >= 0.012}
    lining = random.choice(list(rigid_no_glassy.keys()))
    n = rigid_no_glassy[lining]
    Q = round(random.uniform(5.0, 40.0), 1)
    S = round(random.uniform(0.0005, 0.0018), 4)

    c = round(2 ** (1.0 / 3.0), 4)
    sqS = round(math.sqrt(S), 5)
    y = round((Q * n / (c * sqS)) ** 0.375, 3)
    b = round(2 * y, 3)

    A = round(2 * y ** 2, 3)
    R = round(y / 2, 4)
    Q_check = round((1.0 / n) * A * R ** (2.0 / 3.0) * sqS, 2)
    Fr = (Q / A) / math.sqrt(GRAVITY_M_S2 * y)

    assert 0.9 <= y <= 3.5, f"depth out of bounds: {y}"
    assert Fr <= 0.90, f"best-section flow not subcritical: Fr = {Fr}"
    assert abs(Q_check - Q) / Q <= 0.015, f"back-check failed: {Q_check} vs {Q}"

    question = (
        f"A rectangular channel with {_lining_phrase(lining)} (Manning's "
        f"n = {n}) must carry a design discharge of {Q:.1f} m^3/s on a "
        f"slope of S = {S}. The channel is to be proportioned as the "
        f"best hydraulic (most efficient) rectangular section, for which "
        f"the width is twice the depth. Determine the required flow "
        f"depth y in m."
    )

    solution = (
        f"**Given:**\n"
        f"Design discharge (Q): {Q:.1f} m^3/s; slope (S): {S}; "
        f"n = {n}\n\n"
        f"**Step 1:** State the best hydraulic rectangular section "
        f"properties.\n"
        f"For the most efficient rectangular section, b = 2y, so "
        f"A = 2y^2, P = 4y, and R = A/P = y/2.\n\n"
        f"**Step 2:** Substitute into Manning's equation and collect "
        f"powers of y.\n"
        f"Q = (1/n) * 2y^2 * (y/2)^(2/3) * S^(1/2) "
        f"= 2^(1/3) * y^(8/3) * S^(1/2) / n, with 2^(1/3) = {c:.4f} and "
        f"S^(1/2) = {sqS:.5f}\n\n"
        f"**Step 3:** Solve for the depth.\n"
        f"y = [Q*n / ({c:.4f} * {sqS:.5f})]^(3/8) "
        f"= [{Q:.1f} * {n} / ({c:.4f} * {sqS:.5f})]^(3/8) = {y:.3f} m\n"
        f"(width b = 2y = {b:.3f} m)\n\n"
        f"**Step 4:** Back-check with Manning's equation.\n"
        f"A = 2y^2 = {A:.3f} m^2, R = y/2 = {R:.4f} m\n"
        f"Q = (1/{n}) * {A:.3f} * ({R:.4f})^(2/3) * {sqS:.5f} "
        f"= {Q_check:.2f} m^3/s, matching the design discharge.\n\n"
        f"**Answer:** The required flow depth is {y:.3f} m"
    )

    return question, solution


# Template 24 (Advanced) — Area C1: Uniform Open-Channel Flow
def template_normal_depth_iteration():
    """
    Normal Depth by Trial and Linear Interpolation

    Scenario:
        Manning's equation cannot be inverted in closed form for the
        normal depth of a channel of general shape, so the depth is found
        iteratively. The required section factor is

            A*R^(2/3) = Q*n / S^(1/2) = K   (units m^(8/3))

        The QUESTION prescribes the scheme (so the gold path is
        reproducible): evaluate g(y) = A*R^(2/3) - K at trial depths
        1.000 m and 1.500 m, then update by linear interpolation (secant)
        between successive trials until the depth changes by less than
        0.002 m. Each g-evaluation shows its geometry (A, P, A*R^(2/3)).
        The channel shape (rectangular or trapezoidal) changes the
        geometry functions inside every evaluation.

    Difficulty: Advanced
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 4 (normal
        depth computation; section factor A*R^(2/3)); n per constants
        (HDS-4, on-disk), rigid linings only.
    Physical bounds: target normal depth sampled in [0.8, 2.0] m with the
        discharge DERIVED from it; b in [2.0, 5.0] m; trapezoidal z in
        {1.5, 2.0, 2.5}; slope sampled inside a per-sample subcritical
        window (Fr at the normal depth asserted <= 0.92); convergence
        within 5 interpolation updates (asserted via the tolerance);
        final depth within 0.015 m of the sampled target.

    Returns:
        tuple: (question, solution)
    """
    shape = random.choice(["rectangular", "trapezoidal"])
    lining = random.choice(list(_RIGID_LININGS.keys()))
    n = _RIGID_LININGS[lining]
    b = round(random.uniform(2.0, 5.0), 1)
    z = random.choice([1.5, 2.0, 2.5]) if shape == "trapezoidal" else 0.0
    yn_target = round(random.uniform(0.8, 2.0), 2)

    def geometry(depth):
        if shape == "rectangular":
            A = b * depth
            P = b + 2 * depth
            T = b
        else:
            A = (b + z * depth) * depth
            P = b + 2 * depth * math.sqrt(1 + z ** 2)
            T = b + 2 * z * depth
        return A, P, T

    A_n, P_n, T_n = geometry(yn_target)
    S = _froude_capped_slope(n, A_n, P_n, T_n, 90.0,
                             s_floor=0.0008, s_ceil=0.003)

    def section_factor(depth):
        A, P, _ = geometry(depth)
        return A * (A / P) ** (2.0 / 3.0)

    sqS = round(math.sqrt(S), 5)
    Q = round(section_factor(yn_target) * sqS / n, 2)
    K = round(Q * n / sqS, 3)

    def fmt3(x):
        s = f"{x:.3f}"
        return "0.000" if s == "-0.000" else s

    def evaluate(depth):
        A, P, _ = geometry(depth)
        A_r = round(A, 3)
        P_r = round(P, 3)
        AR = round(section_factor(depth), 3)
        g = round(AR - K, 3)
        return A_r, P_r, AR, g

    eval_lines = []
    y_prev, y_curr = 1.000, 1.500
    A0, P0, AR0, g_prev = evaluate(y_prev)
    eval_lines.append(
        f"At y = {y_prev:.4f} m: A = {A0:.3f} m^2, P = {P0:.3f} m, "
        f"A*R^(2/3) = {AR0:.3f}, so g = {AR0:.3f} - {K:.3f} = "
        f"{fmt3(g_prev)}")
    A1, P1, AR1, g_curr = evaluate(y_curr)
    eval_lines.append(
        f"At y = {y_curr:.4f} m: A = {A1:.3f} m^2, P = {P1:.3f} m, "
        f"A*R^(2/3) = {AR1:.3f}, so g = {AR1:.3f} - {K:.3f} = "
        f"{fmt3(g_curr)}")

    updates = 0
    for _ in range(5):
        y_next = round(
            y_curr - g_curr * (y_curr - y_prev) / (g_curr - g_prev), 4)
        updates += 1
        eval_lines.append(
            f"Update {updates}: y_next = {y_curr:.4f} - ({fmt3(g_curr)}) "
            f"* ({y_curr:.4f} - {y_prev:.4f}) / (({fmt3(g_curr)}) - "
            f"({fmt3(g_prev)})) = {y_next:.4f} m")
        if abs(y_next - y_curr) < 0.002:
            y_curr = y_next
            break
        y_prev, g_prev = y_curr, g_curr
        y_curr = y_next
        A2, P2, AR2, g_curr = evaluate(y_curr)
        eval_lines.append(
            f"At y = {y_curr:.4f} m: A = {A2:.3f} m^2, P = {P2:.3f} m, "
            f"A*R^(2/3) = {AR2:.3f}, so g = {AR2:.3f} - {K:.3f} = "
            f"{fmt3(g_curr)}")
    yn = round(y_curr, 3)
    iter_text = "\n".join(eval_lines)

    A_f, _, _ = geometry(yn)
    Fr = (Q / A_f) / math.sqrt(GRAVITY_M_S2 * (A_f / geometry(yn)[2]))
    assert updates <= 5 and abs(yn - yn_target) <= 0.015, (
        f"iteration failed: {updates} updates, yn = {yn} vs {yn_target}")
    assert abs(section_factor(yn) - K) / K <= 0.02, "residual too large"
    assert Fr <= 0.92, f"normal flow not subcritical: Fr = {Fr}"

    if shape == "rectangular":
        geom_text = f"a rectangular channel of width b = {b:.1f} m"
        geom_note = "A = b*y and P = b + 2y"
    else:
        geom_text = (f"a trapezoidal channel of bottom width b = "
                     f"{b:.1f} m and side slopes {z}H:1V")
        geom_note = (f"A = (b + {z}*y)*y and "
                     f"P = b + 2*y*sqrt(1 + {z}^2)")

    question = (
        f"Uniform flow of Q = {Q:.2f} m^3/s occurs in {geom_text}, with "
        f"{_lining_phrase(lining)} (Manning's n = {n}) on a slope of S = {S}. "
        f"Determine the normal depth of flow in m as follows: form the "
        f"required section factor K = Q*n/S^(1/2), define "
        f"g(y) = A*R^(2/3) - K (A = flow area, P = wetted perimeter, "
        f"R = A/P), evaluate g at trial depths of 1.000 m "
        f"and 1.500 m, and then update the depth by linear interpolation "
        f"(secant) between successive trials until the depth changes by "
        f"less than 0.002 m."
    )

    solution = (
        f"**Given:**\n"
        f"Q = {Q:.2f} m^3/s; {geom_text}; n = {n}; S = {S}\n\n"
        f"**Step 1:** Form the required section factor from Manning's "
        f"equation.\n"
        f"K = A*R^(2/3) = Q*n / S^(1/2) = {Q:.2f} * {n} / {sqS:.5f} "
        f"= {K:.3f} m^(8/3)\n\n"
        f"**Step 2:** Set up the prescribed trial-and-interpolation "
        f"scheme.\n"
        f"With {geom_note}, the section factor A*R^(2/3) increases "
        f"monotonically with depth. Evaluate g(y) = A*R^(2/3) - K at the "
        f"trial depths 1.000 m and 1.500 m, then interpolate linearly "
        f"between successive trials until the depth change is below "
        f"0.002 m.\n\n"
        f"**Step 3:** Iterate.\n"
        f"{iter_text}\n\n"
        f"**Step 4:** State the converged normal depth.\n"
        f"The last change is below the tolerance, so yn = {yn:.3f} m\n\n"
        f"**Answer:** The normal depth is {yn:.3f} m"
    )

    return question, solution
