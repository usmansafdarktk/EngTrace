import math
import random

from pilot.templates.branches.civil_engineering.constants import GRAVITY_M_S2


# Template 25 (Easy) — Area C2: Energy Principles & Rapidly Varied Flow
def template_critical_depth_froude_classification():
    """
    Critical Depth and Froude-Number Classification

    Scenario:
        A rectangular channel carries a known discharge at a known depth.
        The critical depth and the Froude number classify the flow
        regime; the sampled depth lies on either side of critical, so the
        classification conclusion branches:

            q = Q / b;   yc = (q^2 / g)^(1/3)
            Fr = q / (y * sqrt(g*y))
            y > yc  <->  Fr < 1  (subcritical); the converse for
            supercritical.

    Difficulty: Easy
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 2 (specific
        energy; critical depth and Froude number in rectangular
        channels).
    Physical bounds: width b in [2.0, 6.0] m; discharge Q in
        [4, 40] m^3/s; the actual depth is sampled at [1.15, 1.60]*yc
        (subcritical branch) or [0.55, 0.85]*yc (supercritical branch),
        never within 10% of critical; Fr asserted in [0.30, 2.60].

    Returns:
        tuple: (question, solution)
    """
    g = GRAVITY_M_S2
    b = round(random.uniform(2.0, 6.0), 1)
    Q = round(random.uniform(4.0, 40.0), 1)
    regime = random.choice(["subcritical", "supercritical"])

    q = round(Q / b, 3)
    yc = round((q ** 2 / g) ** (1.0 / 3.0), 3)
    if regime == "subcritical":
        y = round(yc * random.uniform(1.15, 1.60), 2)
    else:
        y = round(yc * random.uniform(0.55, 0.85), 2)

    Fr = round(q / (y * math.sqrt(g * y)), 3)

    assert 0.30 <= Fr <= 2.60, f"Froude number out of bounds: {Fr}"
    if regime == "subcritical":
        assert y > yc and Fr < 1.0, f"regime inconsistent: {y}, {yc}, {Fr}"
    else:
        assert y < yc and Fr > 1.0, f"regime inconsistent: {y}, {yc}, {Fr}"

    conclusion = (
        f"Since y = {y:.2f} m > yc = {yc:.3f} m (equivalently Fr < 1), "
        f"the flow is subcritical."
        if regime == "subcritical" else
        f"Since y = {y:.2f} m < yc = {yc:.3f} m (equivalently Fr > 1), "
        f"the flow is supercritical."
    )

    question = (
        f"A rectangular channel {b:.1f} m wide carries a discharge of "
        f"Q = {Q:.1f} m^3/s at a flow depth of {y:.2f} m. Taking "
        f"g = {g} m/s^2, determine the Froude number of the flow. In your "
        f"solution, compute the critical depth and use it to justify "
        f"whether the flow is subcritical or supercritical."
    )

    solution = (
        f"**Given:**\n"
        f"Width (b): {b:.1f} m; discharge (Q): {Q:.1f} m^3/s; depth (y): "
        f"{y:.2f} m; g = {g} m/s^2\n\n"
        f"**Step 1:** Compute the unit discharge.\n"
        f"q = Q / b = {Q:.1f} / {b:.1f} = {q:.3f} m^2/s\n\n"
        f"**Step 2:** Compute the critical depth.\n"
        f"yc = (q^2 / g)^(1/3) = (({q:.3f})^2 / {g})^(1/3) "
        f"= {yc:.3f} m\n\n"
        f"**Step 3:** Compute the Froude number at the actual depth.\n"
        f"Fr = q / (y * sqrt(g*y)) = {q:.3f} / ({y:.2f} * "
        f"sqrt({g} * {y:.2f})) = {Fr:.3f}\n\n"
        f"**Step 4:** Classify the regime.\n"
        f"{conclusion}\n\n"
        f"**Answer:** The Froude number is {Fr:.3f}"
    )

    return question, solution


# Template 26 (Intermediate) — Area C2: Energy Principles & Rapidly Varied Flow
def template_hydraulic_jump_energy_loss():
    """
    Hydraulic Jump: Sequent Depth and Energy Loss

    Scenario:
        A hydraulic jump forms in a rectangular channel. Depending on
        which depth is known — the supercritical upstream depth y1 or the
        subcritical downstream depth y2 — the Belanger relation is
        applied with the corresponding Froude number, an
        equation-selection branch:

            known y1:  Fr1^2 = q^2/(g*y1^3);
                       y2 = (y1/2)*(sqrt(1 + 8*Fr1^2) - 1)
            known y2:  Fr2^2 = q^2/(g*y2^3);
                       y1 = (y2/2)*(sqrt(1 + 8*Fr2^2) - 1)

        The head loss then follows from  dE = (y2 - y1)^3 / (4*y1*y2).

    Difficulty: Intermediate
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 3 (momentum
        principle; hydraulic jump and the Belanger equation).
    Physical bounds: upstream Froude number sampled in [1.8, 5.5]
        (spanning the weak, oscillating, and steady jump classes); y1 sampled inside a per-sample window of
        [0.25, 0.80] m intersected with the head-loss band (dE scales
        linearly with y1 at fixed Fr1); unit discharge derived from
        (Fr1, y1); head loss asserted in [0.05, 2.5] m.

    Returns:
        tuple: (question, solution)
    """
    g = GRAVITY_M_S2
    known = random.choice(["upstream", "downstream"])
    Fr1 = random.uniform(1.8, 5.5)
    # Per-sample y1 window (lesson 1): dE = y1*(r-1)^3/(4r) with
    # r = y2/y1 fixed by Fr1, so bound y1 from the dE band [0.062, 2.38]
    # (buffers sized so rounding can never breach the [0.05, 2.5] assert).
    r_ratio = (math.sqrt(1 + 8 * Fr1 ** 2) - 1) / 2.0
    f_dE = (r_ratio - 1) ** 3 / (4 * r_ratio)
    y1_lo = max(0.25, 0.062 / f_dE)
    y1_hi = min(0.80, 2.38 / f_dE)
    assert y1_lo < y1_hi, f"empty y1 window: {y1_lo}, {y1_hi}"
    y1_true = round(random.uniform(y1_lo, y1_hi), 2)
    q_true = Fr1 * math.sqrt(g * y1_true ** 3)
    y2_true = y1_true / 2 * (math.sqrt(1 + 8 * Fr1 ** 2) - 1)

    q = round(q_true, 3)
    if known == "upstream":
        y1 = y1_true
        Fr1_sq = round(q ** 2 / (g * y1 ** 3), 3)
        y2 = round(y1 / 2 * (math.sqrt(1 + 8 * Fr1_sq) - 1), 4)
        known_text = f"the upstream (supercritical) depth is {y1:.2f} m"
        branch_step = (
            f"**Step 2:** Compute the upstream Froude number.\n"
            f"Fr1^2 = q^2 / (g * y1^3) = ({q:.3f})^2 / ({g} * "
            f"({y1:.2f})^3) = {Fr1_sq:.3f}\n\n"
            f"**Step 3:** Apply the Belanger equation for the sequent "
            f"(downstream) depth.\n"
            f"y2 = (y1/2) * (sqrt(1 + 8*Fr1^2) - 1) "
            f"= ({y1:.2f}/2) * (sqrt(1 + 8 * {Fr1_sq:.3f}) - 1) "
            f"= {y2:.4f} m\n\n"
        )
    else:
        y2 = round(y2_true, 2)
        Fr2_sq = round(q ** 2 / (g * y2 ** 3), 4)
        y1 = round(y2 / 2 * (math.sqrt(1 + 8 * Fr2_sq) - 1), 4)
        known_text = f"the downstream (subcritical) depth is {y2:.2f} m"
        branch_step = (
            f"**Step 2:** Compute the downstream Froude number.\n"
            f"Fr2^2 = q^2 / (g * y2^3) = ({q:.3f})^2 / ({g} * "
            f"({y2:.2f})^3) = {Fr2_sq:.4f}\n\n"
            f"**Step 3:** Apply the Belanger equation in its conjugate "
            f"form for the upstream depth.\n"
            f"y1 = (y2/2) * (sqrt(1 + 8*Fr2^2) - 1) "
            f"= ({y2:.2f}/2) * (sqrt(1 + 8 * {Fr2_sq:.4f}) - 1) "
            f"= {y1:.4f} m\n\n"
        )

    dE = round((y2 - y1) ** 3 / (4 * y1 * y2), 3)
    # Display precisions follow which depth was given (2 dp) vs computed
    # (4 dp — cubed-difference amplification demands the extra digits).
    y1_disp = f"{y1:.2f}" if known == "upstream" else f"{y1:.4f}"
    y2_disp = f"{y2:.4f}" if known == "upstream" else f"{y2:.2f}"

    assert y2 > y1 > 0, f"depth ordering invalid: {y1}, {y2}"
    assert 0.05 <= dE <= 2.5, f"energy loss out of bounds: {dE}"

    question = (
        f"A hydraulic jump forms in a rectangular channel carrying a "
        f"unit discharge of q = {q:.3f} m^2/s ({'per meter of width'}). "
        f"Measurements show that {known_text}. Taking g = {g} m/s^2, "
        f"use the momentum (Belanger) relation to find the sequent "
        f"depth, and then determine the head loss across the jump in m."
    )

    solution = (
        f"**Given:**\n"
        f"Unit discharge (q): {q:.3f} m^2/s; {known_text}; g = {g} "
        f"m/s^2\n\n"
        f"**Step 1:** Identify the known side of the jump.\n"
        f"The jump connects a supercritical upstream depth y1 to its "
        f"subcritical sequent depth y2; here {known_text}, so the "
        f"Belanger relation is applied with the Froude number on the "
        f"known side.\n\n"
        f"{branch_step}"
        f"**Step 4:** Compute the head loss across the jump.\n"
        f"dE = (y2 - y1)^3 / (4*y1*y2) = ({y2_disp} - {y1_disp})^3"
        f" / (4 * {y1_disp} * {y2_disp}) = {dE:.3f} m\n\n"
        f"**Answer:** The head loss across the jump is {dE:.3f} m"
    )

    return question, solution


# Template 27 (Intermediate) — Area C2: Energy Principles & Rapidly Varied Flow
def template_max_hump_height_no_choking():
    """
    Maximum Hump Height Without Choking

    Scenario:
        Subcritical flow approaches a smooth bed hump. The flow chokes
        when the specific energy over the hump falls to the critical
        minimum, so the maximum hump height that passes the flow without
        choking is the available energy margin:

            E1 = y1 + q^2/(2*g*y1^2)
            yc = (q^2/g)^(1/3);  Ec = 1.5*yc
            dz_max = E1 - Ec

    Difficulty: Intermediate
    Grounding: Sturm, Open Channel Hydraulics, 1st ed., Ch. 2 (specific
        energy; transitions and choking).
    Physical bounds: critical depth sampled in [0.5, 1.4] m with the
        approach depth at [1.30, 2.00]*yc (firmly subcritical); q derived
        from yc; dz_max asserted in [0.05, 1.2] m.

    Returns:
        tuple: (question, solution)
    """
    g = GRAVITY_M_S2
    yc_target = round(random.uniform(0.5, 1.4), 3)
    q = round(math.sqrt(g * yc_target ** 3), 3)

    # Per-sample multiplier floor (R1, cycle 1: yc near 0.5 with the
    # multiplier near 1.30 gave dz below the 0.05 floor). dz = yc*f(m)
    # with f(m) = m + 1/(2m^2) - 1.5 increasing for m > 1; bisect for the
    # m that gives dz = 0.058, then sample above it (cap 2.00).
    def _f(m):
        return m + 1.0 / (2.0 * m ** 2) - 1.5

    target = 0.058 / yc_target
    m_lo, m_hi_b = 1.30, 2.00
    if _f(m_lo) < target:
        lo, hi = m_lo, m_hi_b
        for _ in range(40):
            mid = (lo + hi) / 2
            if _f(mid) < target:
                lo = mid
            else:
                hi = mid
        m_lo = hi
    y1 = round(yc_target * random.uniform(m_lo, 2.00), 2)

    # Energies carried at 4 dp: dz is a small difference of near-equal
    # energies, and 3-dp rounding of E1/Ec flipped the final digit
    # (R2, cycle 1).
    E1 = round(y1 + q ** 2 / (2 * g * y1 ** 2), 4)
    yc = round((q ** 2 / g) ** (1.0 / 3.0), 3)
    Ec = round(1.5 * yc, 4)
    dz_max = round(E1 - Ec, 3)

    assert y1 > yc, f"approach flow not subcritical: {y1} vs {yc}"
    assert 0.05 <= dz_max <= 1.2, f"hump margin out of bounds: {dz_max}"

    question = (
        f"Water flows at a unit discharge of q = {q:.3f} m^2/s in a wide "
        f"rectangular channel at an approach depth of y1 = {y1:.2f} m "
        f"(subcritical). The bed rises gradually over a smooth hump. "
        f"Taking g = {g} m/s^2 and neglecting energy losses, determine "
        f"the maximum hump height, in m, for which the flow passes "
        f"without choking (i.e., without forcing critical conditions on "
        f"the hump)."
    )

    solution = (
        f"**Given:**\n"
        f"Unit discharge (q): {q:.3f} m^2/s; approach depth (y1): "
        f"{y1:.2f} m; g = {g} m/s^2\n\n"
        f"**Step 1:** Compute the approach specific energy.\n"
        f"E1 = y1 + q^2/(2*g*y1^2) = {y1:.2f} + ({q:.3f})^2 / (2 * {g} "
        f"* ({y1:.2f})^2) = {E1:.4f} m\n\n"
        f"**Step 2:** Compute the critical depth.\n"
        f"yc = (q^2/g)^(1/3) = (({q:.3f})^2 / {g})^(1/3) = {yc:.3f} m\n\n"
        f"**Step 3:** Compute the minimum specific energy at critical "
        f"conditions.\n"
        f"Ec = 1.5 * yc = 1.5 * {yc:.3f} = {Ec:.4f} m\n\n"
        f"**Step 4:** The flow chokes when the energy over the hump "
        f"drops to Ec, so the maximum bed rise equals the energy "
        f"margin.\n"
        f"dz_max = E1 - Ec = {E1:.4f} - {Ec:.4f} = {dz_max:.3f} m\n\n"
        f"**Answer:** The maximum hump height without choking is "
        f"{dz_max:.3f} m"
    )

    return question, solution
