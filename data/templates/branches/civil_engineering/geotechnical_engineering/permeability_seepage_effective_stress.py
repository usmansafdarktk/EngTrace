import math
import random

from pilot.templates.branches.civil_engineering.constants import (
    PERMEABILITY_RANGES_CM_S,
    SPECIFIC_GRAVITY_RANGES,
    UNIT_WEIGHT_WATER_KN_M3,
)


# Template 4 (Easy) — Area B2: Permeability, Seepage & Effective Stress
def template_constant_head_permeability():
    """
    Coefficient of Permeability from a Constant-Head Test

    Scenario:
        A constant-head permeability test on a sand specimen collects a
        measured volume of water over a known time under a constant head
        difference. Darcy's law applied to the test gives

            k = V * L / (A * h * t)

        with A the specimen cross-sectional area.

    Difficulty: Easy
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Ch. 7 (laboratory determination of k — constant-head test). Sampled
        k kept inside the coarse-sand band of Das Table 7.1 as transcribed
        in constants.PERMEABILITY_RANGES_CM_S.
    Physical bounds: specimen D in [5.0, 10.0] cm, L in [10.0, 20.0] cm,
        head h in [20, 60] cm capped at 3L (gradient i <= 3), time t in
        [60, 300] s; underlying k sampled with per-sample feasibility bounds
        inside [0.02, 0.09] cm/s (coarse sand); recomputed k asserted inside
        the Das coarse-sand band [0.01, 1.0] cm/s; collected volume in
        [50, 5000] whole cm^3.

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize the test setup; derive the collected volume from a
    # target k so every instance is a physically consistent test record.
    D = round(random.uniform(5.0, 10.0), 1)      # specimen diameter, cm
    L = round(random.uniform(10.0, 20.0), 1)     # specimen length, cm
    # Head capped at 3L so the gradient stays <= 3 (R1, cycle 1: textbook
    # practice; unbounded h/L reached i = 6, outside Darcy lab validity).
    h = random.randint(20, min(60, int(3 * L)))  # constant head, cm
    t = random.randint(6, 30) * 10               # collection time, s

    # Per-sample feasibility bounds on k (AUTHOR_NOTES lesson 1): the
    # derived volume V = k*A*h*t/L must land inside [50, 5000] cm^3 for
    # EVERY joint draw of the test geometry, so k's sampling window is
    # computed from the geometry (then intersected with the coarse-sand
    # band). Feasibility holds across the whole geometry space: the
    # tightest window is [0.02, 0.0347] at the largest-A/h/t, smallest-L
    # corner.
    A_true = math.pi * D ** 2 / 4.0
    k_lo_s = max(0.02, 1.02 * 50.0 * L / (A_true * h * t))
    k_hi_s = min(0.09, 0.98 * 5000.0 * L / (A_true * h * t))
    k_true = random.uniform(k_lo_s, k_hi_s)      # coarse sand, cm/s

    V = round(k_true * A_true * h * t / L, 0)    # presented volume, whole cm^3

    # 2. Core computation — round-then-recompute from presented values.
    A = round(math.pi * D ** 2 / 4.0, 2)
    k = round(V * L / (A * h * t), 4)

    k_lo, k_hi = PERMEABILITY_RANGES_CM_S["coarse sand"]
    assert k_lo <= k <= k_hi, f"k outside Das coarse-sand band: {k}"
    assert 50.0 <= V <= 5000.0, f"collected volume implausible: {V}"

    # 3. Serialize.
    question = (
        f"A constant-head permeability test is run on a specimen of coarse "
        f"sand. The cylindrical specimen has a diameter of {D:.1f} cm and a "
        f"length of {L:.1f} cm, and the head difference across the specimen "
        f"is held constant at {h} cm. Over a period of {t} seconds, "
        f"{V:.0f} cm^3 of water is collected at the outlet. Determine the "
        f"coefficient of permeability of the sand in cm/s."
    )

    solution = (
        f"**Given:**\n"
        f"Specimen diameter (D): {D:.1f} cm\n"
        f"Specimen length (L): {L:.1f} cm\n"
        f"Constant head difference (h): {h} cm\n"
        f"Collection time (t): {t} s\n"
        f"Collected volume (V): {V:.0f} cm^3\n\n"
        f"**Step 1:** Compute the specimen cross-sectional area.\n"
        f"A = pi * D^2 / 4 = pi * ({D:.1f})^2 / 4 = {A:.2f} cm^2\n\n"
        f"**Step 2:** Apply Darcy's law for the constant-head test.\n"
        f"For steady flow, V = k * i * A * t with i = h / L, which gives\n"
        f"k = V * L / (A * h * t)\n\n"
        f"**Step 3:** Substitute and evaluate.\n"
        f"k = ({V:.0f} * {L:.1f}) / ({A:.2f} * {h} * {t}) = {k:.4f} cm/s\n\n"
        f"**Answer:** The coefficient of permeability is {k:.4f} cm/s"
    )

    return question, solution


# Template 5 (Intermediate) — Area B2: Permeability, Seepage & Effective Stress
def template_effective_stress_profile():
    """
    Effective Stress Below a Water Table

    Scenario:
        A sand deposit has its water table at depth z1. Above the water
        table the sand is moist (dry unit weight and moisture content
        given); below it the sand is saturated (void ratio and specific
        gravity given). The chain is:

            gamma_moist = gamma_d * (1 + w)
            gamma_sat   = (Gs + e) * gamma_w / (1 + e)
            sigma  = gamma_moist * z1 + gamma_sat * z2
            u      = gamma_w * z2
            sigma' = sigma - u

    Difficulty: Intermediate
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Ch. 9 (in situ stresses; effective stress without seepage), with
        phase relationships from Ch. 3.
    Physical bounds: above the water table gamma_d in [14.5, 16.5] kN/m^3
        and w in [6%, 14%]; below it Gs per sand and e in [0.45,
        min(0.75, e_implied_above - 0.02)] so the deposit never loosens
        with depth and gamma_sat always exceeds gamma_moist;
        z1 in [2.0, 5.0] m, z2 in [3.0, 8.0] m; effective stress asserted
        in [40, 220] kPa.

    Returns:
        tuple: (question, solution)
    """
    gamma_w = UNIT_WEIGHT_WATER_KN_M3

    # 1. Parameterize the two zones independently (ranges chosen so the
    # saturated zone is always at least as heavy as the moist zone).
    gamma_d = round(random.uniform(14.5, 16.5), 2)   # above WT, kN/m^3
    w_pct = round(random.uniform(6.0, 14.0), 1)      # above WT, %
    gs_lo, gs_hi = SPECIFIC_GRAVITY_RANGES["sand"]
    Gs = round(random.uniform(gs_lo, gs_hi), 2)      # deposit-wide
    # Couple the zones (R1+R3, cycle 1): the saturated zone must be at
    # least as dense as the moist zone above it, so its void ratio is capped
    # by the void ratio implied by the moist zone's gamma_d.
    e_above = Gs * gamma_w / gamma_d - 1
    e = round(random.uniform(0.45, min(0.75, e_above - 0.02)), 2)
    z1 = round(random.uniform(2.0, 5.0), 1)          # depth to WT, m
    z2 = round(random.uniform(3.0, 8.0), 1)          # WT to point A, m

    # 2. Core computation — round-then-recompute at every step.
    w = round(w_pct / 100.0, 3)
    gamma_moist = round(gamma_d * (1 + w), 2)
    gamma_sat = round((Gs + e) * gamma_w / (1 + e), 2)
    sigma = round(gamma_moist * z1 + gamma_sat * z2, 1)
    u = round(gamma_w * z2, 1)
    sigma_eff = round(sigma - u, 1)

    assert gamma_moist < gamma_sat, (
        f"moist zone heavier than saturated zone: {gamma_moist}, {gamma_sat}")
    assert 40.0 <= sigma_eff <= 220.0, (
        f"effective stress out of bounds: {sigma_eff}")

    # 3. Serialize.
    question = (
        f"At a site, a deep deposit of sand has its groundwater table at a "
        f"depth of {z1:.1f} m below the ground surface. Above the water "
        f"table the sand is moist, with a dry unit weight of "
        f"{gamma_d:.2f} kN/m^3 and a moisture content of {w_pct:.1f}%. "
        f"Below the water table the sand is saturated, with a void ratio of "
        f"{e:.2f} and a specific gravity of solids of {Gs:.2f}. Taking the "
        f"unit weight of water as {gamma_w:.2f} kN/m^3, determine the "
        f"vertical effective stress in kPa at point A, located "
        f"{z2:.1f} m below the water table."
    )

    solution = (
        f"**Given:**\n"
        f"Depth to water table (z1): {z1:.1f} m\n"
        f"Moist zone: gamma_d = {gamma_d:.2f} kN/m^3, w = {w_pct:.1f}% "
        f"= {w:.3f}\n"
        f"Saturated zone: e = {e:.2f}, Gs = {Gs:.2f}\n"
        f"Depth of A below water table (z2): {z2:.1f} m\n"
        f"Unit weight of water (gamma_w): {gamma_w:.2f} kN/m^3\n\n"
        f"**Step 1:** Compute the moist unit weight above the water table.\n"
        f"gamma_moist = gamma_d * (1 + w) = {gamma_d:.2f} * (1 + {w:.3f}) "
        f"= {gamma_moist:.2f} kN/m^3\n\n"
        f"**Step 2:** Compute the saturated unit weight below the water "
        f"table.\n"
        f"gamma_sat = (Gs + e) * gamma_w / (1 + e) "
        f"= ({Gs:.2f} + {e:.2f}) * {gamma_w:.2f} / (1 + {e:.2f}) "
        f"= {gamma_sat:.2f} kN/m^3\n\n"
        f"**Step 3:** Compute the total vertical stress at A.\n"
        f"sigma = gamma_moist * z1 + gamma_sat * z2 "
        f"= {gamma_moist:.2f} * {z1:.1f} + {gamma_sat:.2f} * {z2:.1f} "
        f"= {sigma:.1f} kPa\n\n"
        f"**Step 4:** Compute the pore water pressure at A.\n"
        f"u = gamma_w * z2 = {gamma_w:.2f} * {z2:.1f} = {u:.1f} kPa\n\n"
        f"**Step 5:** Apply the effective stress principle.\n"
        f"sigma' = sigma - u = {sigma:.1f} - {u:.1f} = {sigma_eff:.1f} kPa\n\n"
        f"**Answer:** The vertical effective stress at point A is "
        f"{sigma_eff:.1f} kPa"
    )

    return question, solution


# Template 6 (Intermediate) — Area B2: Permeability, Seepage & Effective Stress
def template_upward_seepage_quick_condition():
    """
    Factor of Safety Against the Quick Condition Under Upward Seepage

    Scenario:
        Upward seepage flows through a sand layer of thickness L under a
        measured excess head h across the layer. The exit hydraulic
        gradient is compared with the critical gradient at which effective
        stress vanishes:

            i     = h / L
            gamma' = (Gs - 1) * gamma_w / (1 + e)
            i_cr  = gamma' / gamma_w
            FS    = i_cr / i

    Difficulty: Intermediate
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Section 9.3 (Stresses in Saturated Soil with Upward Seepage;
        critical hydraulic gradient and the quick/boiling condition).
    Physical bounds: Gs per sand; e in [0.45, 0.85] so i_cr in ~[0.89,
        1.15]; target FS sampled in [1.25, 3.0] (upward seepage present but
        stable); layer thickness L in [1.5, 4.0] m; derived excess head h
        in [0.45, 3.69] m; i carried at 4 decimals (precision sized to
        downstream division).

    Returns:
        tuple: (question, solution)
    """
    gamma_w = UNIT_WEIGHT_WATER_KN_M3

    # 1. Parameterize: sample the soil state and a target FS, derive the
    # excess head so every instance is jointly consistent.
    gs_lo, gs_hi = SPECIFIC_GRAVITY_RANGES["sand"]
    Gs = round(random.uniform(gs_lo, gs_hi), 2)
    e = round(random.uniform(0.45, 0.85), 2)
    L = round(random.uniform(1.5, 4.0), 1)
    FS_true = random.uniform(1.25, 3.0)

    i_cr_true = (Gs - 1) / (1 + e)
    h = round(i_cr_true / FS_true * L, 2)        # presented excess head, m

    # 2. Core computation — round-then-recompute at every step; i and i_cr
    # carry 4 decimals (they feed a division; see AUTHOR_NOTES lesson 5).
    i = round(h / L, 4)
    gamma_sub = round((Gs - 1) * gamma_w / (1 + e), 2)
    i_cr = round(gamma_sub / gamma_w, 4)
    FS = round(i_cr / i, 3)

    assert 0.85 <= i_cr <= 1.16, f"critical gradient out of bounds: {i_cr}"
    assert 1.15 <= FS <= 3.15, f"factor of safety out of bounds: {FS}"
    # Ceiling dominates the analytic max of the derived head:
    # i_cr_max/FS_min * L_max = 1.1517/1.25 * 4.0 = 3.69 (R1, cycle 1).
    assert 0.30 <= h <= 3.70, f"excess head implausible: {h}"

    # 3. Serialize.
    question = (
        f"An excavation is underlain by a {L:.1f} m thick layer of sand "
        f"through which water seeps vertically upward. Piezometer readings "
        f"show that the head loss across the layer is {h:.2f} m. The sand "
        f"has a void ratio of {e:.2f} and a specific gravity of solids of "
        f"{Gs:.2f}, and the unit weight of water is {gamma_w:.2f} kN/m^3. "
        f"Determine the factor of safety against the quick (boiling) "
        f"condition."
    )

    solution = (
        f"**Given:**\n"
        f"Sand layer thickness (L): {L:.1f} m\n"
        f"Head loss across the layer (h): {h:.2f} m\n"
        f"Void ratio (e): {e:.2f}\n"
        f"Specific gravity of solids (Gs): {Gs:.2f}\n"
        f"Unit weight of water (gamma_w): {gamma_w:.2f} kN/m^3\n\n"
        f"**Step 1:** Compute the upward hydraulic gradient through the "
        f"layer.\n"
        f"i = h / L = {h:.2f} / {L:.1f} = {i:.4f}\n\n"
        f"**Step 2:** Compute the submerged (effective) unit weight of the "
        f"sand.\n"
        f"gamma' = (Gs - 1) * gamma_w / (1 + e) "
        f"= ({Gs:.2f} - 1) * {gamma_w:.2f} / (1 + {e:.2f}) "
        f"= {gamma_sub:.2f} kN/m^3\n\n"
        f"**Step 3:** Compute the critical hydraulic gradient.\n"
        f"The quick condition occurs when the upward seepage force cancels "
        f"the submerged weight:\n"
        f"i_cr = gamma' / gamma_w = {gamma_sub:.2f} / {gamma_w:.2f} "
        f"= {i_cr:.4f}\n\n"
        f"**Step 4:** Compute the factor of safety.\n"
        f"FS = i_cr / i = {i_cr:.4f} / {i:.4f} = {FS:.3f}\n\n"
        f"**Answer:** The factor of safety against the quick condition is "
        f"{FS:.3f}"
    )

    return question, solution
