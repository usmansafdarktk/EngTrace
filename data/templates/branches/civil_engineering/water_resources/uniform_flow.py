import math
import random
from fractions import Fraction

from data.templates.branches.civil_engineering.constants import (
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


# --- Display-tie screening (D-016), used by template_normal_depth_iteration --
#
# A printed line closes for a decimal reader and for a binary reader alike only
# when the exact value of its expression is NOT half-way between two displayable
# values. At such a tie |evaluated - printed| == tol exactly and float
# representation error alone decides which side it lands on, in either rounding
# direction -- so the instance has no defensible gold reading and is REMOVED
# rather than resolved (D-016). D-037's cheaper escape, lengthening the display
# until the quantity is exact, does not apply to either quantity screened here:
# both are quotients by a value that is not a power of ten.

def _exact(x, spec):
    """The exact rational value of `x` as the trace prints it under `spec`."""
    return Fraction(format(x, spec))


def _is_display_tie(value, places):
    """Does an exact rational land exactly on a half-way display boundary?"""
    return (value * 10 ** places) % 1 == Fraction(1, 2)


# Iteration budget and convergence tolerance for the normal-depth secant
# scheme. Both appear in the QUESTION text, so they are part of the item's
# statement, not tuning knobs: _T24_TOL is what the solver is told to iterate
# to, and _T24_MAX_UPDATES is a guard whose breach raises rather than
# truncating the trace.
_T24_TOL = 0.002
_T24_MAX_UPDATES = 5


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

    Trace shape (Phase 3, D3.1 — schema route):
        The number of interpolation updates is DATA-DEPENDENT: measured
        over 4,000 seeds it is {1: 29, 2: 312, 3: 2345, 4: 1261, 5: 53},
        so 32.9% of instances need MORE than three updates, because the
        question's termination condition is a tolerance on the depth
        change, not a fixed pass count. The count is therefore INCIDENTAL
        to the answer — a solver reaching the same converged depth in four
        updates instead of three is correct. This template
        builds an `iteration` trace node (the local `trace_nodes`,
        specified in docs/re-implementation-sep/phase3_node_types.md) and
        RENDERS the printed trace from it, so the structured trace and the
        prose cannot disagree.

    Difficulty: Advanced
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 4 (normal
        depth computation; section factor A*R^(2/3)); n per constants
        (HDS-4, on-disk), rigid linings only.
    Physical bounds: target normal depth sampled in [0.8, 2.0] m with the
        discharge DERIVED from it; b in [2.0, 5.0] m; trapezoidal z in
        {1.5, 2.0, 2.5}; slope sampled inside a per-sample subcritical
        window (Fr at the normal depth asserted <= 0.92); convergence
        within 5 interpolation updates — enforced by an explicit raise,
        NOT an assert, because Step 4 tells the reader the tolerance was
        met and `python -O` strips asserts (Phase 2 Reviewer C, C-6);
        final depth within 0.015 m of the sampled target. Screens
        (bounded resample loop, cap 200): an instance is rejected when
        Step 1's K expression or any update expression lands EXACTLY on a
        half-way display boundary, which is an ill-posed instance rather
        than a rounding-convention choice (D-016). Measured rejection
        rate 2.14% over 20,000 seeds (K line 1.14%, update lines 1.02%).

    Returns:
        tuple: (question, solution)
    """
    for _attempt in range(200):
        shape = random.choice(["rectangular", "trapezoidal"])
        lining = random.choice(list(_RIGID_LININGS.keys()))
        n = _RIGID_LININGS[lining]
        b = round(random.uniform(2.0, 5.0), 1)
        z = random.choice([1.5, 2.0, 2.5]) if shape == "trapezoidal" else 0.0
        yn_target = round(random.uniform(0.8, 2.0), 2)

        def geometry(depth, _b=b, _z=z, _shape=shape):
            if _shape == "rectangular":
                A = _b * depth
                P = _b + 2 * depth
                T = _b
            else:
                A = (_b + _z * depth) * depth
                P = _b + 2 * depth * math.sqrt(1 + _z ** 2)
                T = _b + 2 * _z * depth
            return A, P, T

        A_n, P_n, T_n = geometry(yn_target)
        S = _froude_capped_slope(n, A_n, P_n, T_n, 90.0,
                                 s_floor=0.0008, s_ceil=0.003)

        def section_factor(depth, _geom=geometry):
            A, P, _T = _geom(depth)
            return A * (A / P) ** (2.0 / 3.0)

        sqS = round(math.sqrt(S), 5)
        Q = round(section_factor(yn_target) * sqS / n, 2)
        K = round(Q * n / sqS, 3)

        # Step 1's printed line, read exactly in decimal. A tie here is an
        # ill-posed instance, not a rounding-convention choice (D-016).
        # D-037's "lengthen the display" escape was tested first and does not
        # apply: it needs the quantity exactly representable at the longer
        # display for EVERY instance, and K = Q*n/sqrt(S) is a quotient by a
        # 5-dp square root that terminates at 3 dp on only 10.9% of instances,
        # at 4 dp on 11.3% and at 6 dp on 12.9% (3,000 seeds). Rounding
        # persists at any length, so the ties would relocate, not vanish.
        if _is_display_tie(
                _exact(Q, ".2f") * Fraction(str(n)) / _exact(sqS, ".5f"), 3):
            continue

        def evaluate(depth, _K=K, _geom=geometry):
            """One g-evaluation frame: the `evaluation` sub-node of D3.3."""
            A, P, _T = _geom(depth)
            AR = round(A * (A / P) ** (2.0 / 3.0), 3)
            g = round(AR - _K, 3)
            if g == 0.0:
                # Normalise -0.0. The removed fmt3() helper rewrote the
                # STRING "-0.000" to "0.000", which let the printed operand
                # differ in sign from the stored one (a P2 violation).
                # Normalising the stored value keeps the two identical by
                # construction instead of patching the display.
                g = 0.0
            return {"y": round(depth, 4), "A": round(A, 3),
                    "P": round(P, 3), "AR": AR, "g": g}

        # ---- build the `iteration` node; the prose is rendered from it ----
        preamble = [evaluate(1.0000), evaluate(1.5000)]
        elements = []
        y_prev, y_curr = 1.0000, 1.5000
        g_prev, g_curr = preamble[0]["g"], preamble[1]["g"]
        converged = False
        display_tie = False
        for k in range(1, _T24_MAX_UPDATES + 1):
            denom = g_curr - g_prev
            if denom == 0.0:
                # The secant step is undefined. Not observed in 20,000 seeds
                # (smallest |g_k - g_(k-1)| seen is 0.004, four steps of the
                # 3-dp residual grid), so this raises rather than resampling:
                # if the sampling box ever moves it must not pass quietly.
                raise ZeroDivisionError(
                    f"secant denominator vanished at update {k}: "
                    f"g_prev = g_curr = {g_curr:.3f}")
            y_c, y_p = _exact(y_curr, ".4f"), _exact(y_prev, ".4f")
            g_c, g_p = _exact(g_curr, ".3f"), _exact(g_prev, ".3f")
            if _is_display_tie(y_c - g_c * (y_c - y_p) / (g_c - g_p), 4):
                display_tie = True
                break
            y_next = round(y_curr - g_curr * (y_curr - y_prev) / denom, 4)
            change = round(abs(y_next - y_curr), 4)
            element = {"k": k,
                       "y_prev": y_prev, "g_prev": g_prev,
                       "y_curr": y_curr, "g_curr": g_curr,
                       "y_next": y_next, "change": change,
                       "converged": change < _T24_TOL,
                       "evaluation": None}
            elements.append(element)
            if element["converged"]:
                y_curr = y_next
                converged = True
                break
            y_prev, g_prev = y_curr, g_curr
            y_curr = y_next
            element["evaluation"] = evaluate(y_curr)
            g_curr = element["evaluation"]["g"]
        if display_tie:
            continue
        if not converged:
            # NOT an assert: Step 4 tells the reader the tolerance was met,
            # and `python -O` strips asserts (Phase 2 Reviewer C, C-6).
            raise RuntimeError(
                f"secant iteration did not converge within "
                f"{_T24_MAX_UPDATES} updates: last change "
                f"{elements[-1]['change']:.4f} m >= {_T24_TOL} m")
        break
    else:
        raise AssertionError("resample loop exhausted")

    # The `iteration` node. Every name an external verifier needs is DECLARED
    # here rather than inferred: `roles` binds this template's own symbols to
    # the type-level roles of D3.3 §4, `update_relation` is the recurrence as
    # data over those roles, and `symbol_precision` gives the display precision
    # the comparator's tolerance is taken from. Phase 3 Reviewer D showed that
    # without these a verifier can only be written by hardcoding this
    # template's symbol names and its secant formula (finding F3).
    trace_nodes = {
        "schema_version": "1.2",
        "node_type": "iteration",
        "node_id": "t24_secant_normal_depth",
        "cardinality": "incidental",
        "roles": {"index": "k",
                  "iterate_prev": "y_prev", "iterate_curr": "y_curr",
                  "iterate_next": "y_next",
                  "residual_prev": "g_prev", "residual_curr": "g_curr",
                  "change": "change", "converged": "converged",
                  "evaluation": "evaluation"},
        "update_relation": ("iterate_curr - residual_curr "
                            "* (iterate_curr - iterate_prev) "
                            "/ (residual_curr - residual_prev)"),
        "termination": {"kind": "convergence",
                        "quantity_role": "change",
                        "comparison": "lt",
                        "tolerance": _T24_TOL,
                        "predicate_prose": "abs(y_next - y_curr) < 0.002",
                        "max_elements": _T24_MAX_UPDATES,
                        "satisfied": True},
        "rounding": "decimal-half-up",
        "symbol_precision": {"y_prev": 4, "y_curr": 4, "y_next": 4,
                             "change": 4, "g_prev": 3, "g_curr": 3,
                             "y": 4, "A": 3, "P": 3, "AR": 3, "g": 3},
        "preamble_binding": {"iterate_prev": 0, "residual_prev": 0,
                             "iterate_curr": 1, "residual_curr": 1},
        "preamble_symbols": ["y", "A", "P", "AR", "g"],
        "evaluation_symbols": ["y", "A", "P", "AR", "g"],
        "evaluation_roles": {"iterate": "y", "residual": "g"},
        "element_symbols": ["k", "y_prev", "g_prev", "y_curr", "g_curr",
                            "y_next", "change", "converged", "evaluation"],
        "carry": {"y_prev": "y_curr", "g_prev": "g_curr",
                  "y_curr": "y_next", "g_curr": "evaluation.g"},
        "carry_notes": {},
        "preamble": preamble,
        "elements": elements,
        "result": {"symbol": "yn", "value": round(y_curr, 3),
                   "unit": "m", "dp": 3, "from": "last.iterate_next"},
    }

    yn = trace_nodes["result"]["value"]
    updates = len(elements)
    A_f, _P_f, T_f = geometry(yn)
    Fr = (Q / A_f) / math.sqrt(GRAVITY_M_S2 * (A_f / T_f))
    assert updates <= _T24_MAX_UPDATES and abs(yn - yn_target) <= 0.015, (
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

    # ---- render the printed trace FROM the node, never alongside it ----
    def _eval_line(frame, _K=K):
        return (f"At y = {frame['y']:.4f} m: A = {frame['A']:.3f} m^2, "
                f"P = {frame['P']:.3f} m, A*R^(2/3) = {frame['AR']:.3f}, "
                f"so g = {frame['AR']:.3f} - {_K:.3f} = {frame['g']:.3f}")

    eval_lines = [_eval_line(frame) for frame in preamble]
    for element in elements:
        y_cur_s, y_prv_s = element["y_curr"], element["y_prev"]
        g_cur_s, g_prv_s = element["g_curr"], element["g_prev"]
        eval_lines.append(
            f"Update {element['k']}: y_next = {y_cur_s:.4f} - "
            f"({g_cur_s:.3f}) * ({y_cur_s:.4f} - {y_prv_s:.4f}) / "
            f"(({g_cur_s:.3f}) - ({g_prv_s:.3f})) "
            f"= {element['y_next']:.4f} m")
        if element["evaluation"] is not None:
            eval_lines.append(_eval_line(element["evaluation"]))
    iter_text = "\n".join(eval_lines)
    last_change = elements[-1]["change"]

    question = (
        f"Uniform flow of Q = {Q:.2f} m^3/s occurs in {geom_text}, with "
        f"{_lining_phrase(lining)} (Manning's n = {n}) on a slope of S = {S}. "
        f"Determine the normal depth of flow in m as follows: form the "
        f"required section factor K = Q*n/S^(1/2), define "
        f"g(y) = A*R^(2/3) - K (A = flow area, P = wetted perimeter, "
        f"R = A/P), evaluate g at trial depths of 1.000 m "
        f"and 1.500 m, and then update the depth by linear interpolation "
        f"(secant) between successive trials until the depth changes by "
        f"less than {_T24_TOL} m."
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
        f"{_T24_TOL} m.\n\n"
        f"**Step 3:** Iterate. The scheme runs until the tolerance is met, "
        f"so the number of updates is not fixed in advance; here it takes "
        f"{updates}.\n"
        f"{iter_text}\n\n"
        f"**Step 4:** State the converged normal depth.\n"
        f"The last change is {last_change:.4f} m, below the tolerance of "
        f"{_T24_TOL} m, so yn = {yn:.3f} m\n\n"
        f"**Answer:** The normal depth is {yn:.3f} m"
    )

    return question, solution
