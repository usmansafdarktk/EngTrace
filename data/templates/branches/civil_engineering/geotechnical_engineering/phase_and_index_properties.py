import random

from pilot.templates.branches.civil_engineering.constants import (
    SPECIFIC_GRAVITY_RANGES,
    UNIT_WEIGHT_WATER_KN_M3,
)


# Template 1 (Easy) — Area B1: Phase Relationships & Index Properties
def template_phase_relations_degree_of_saturation():
    """
    Degree of Saturation from Moist Unit Weight

    Scenario:
        A moist soil sample has a known total (moist) unit weight, moisture
        content, and specific gravity of solids. The fundamental
        weight-volume relationships give, in sequence:

            gamma_d = gamma / (1 + w)          (dry unit weight)
            e = (Gs * gamma_w / gamma_d) - 1   (void ratio)
            S = w * Gs / e                     (degree of saturation)

        The gold trace applies round-then-recompute at every step: each
        printed equation evaluates exactly with its displayed operands.

    Difficulty: Easy
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Sections 3.2-3.3 (weight-volume relationships; Eqs. 3.9, 3.11 and
        the S*e = w*Gs relation). Gs per SPECIFIC_GRAVITY_RANGES; void-ratio
        ranges conditioned on soil type, anchored to Das Table 3.1
        natural-state values.
    Physical bounds: 2.60 <= Gs <= 2.80; e sampled per soil type (sand
        0.45-0.85, silt 0.50-0.90, inorganic clay 0.55-0.95); S sampled with
        per-sample feasibility bounds so presented w stays in [5%, 33%];
        recomputed S in [15%, 100%]; presented moist unit weight in
        [14, 23] kN/m^3.

    Returns:
        tuple: (question, solution)
    """
    gamma_w = UNIT_WEIGHT_WATER_KN_M3

    # 1. Parameterize: sample the independent phase variables, derive the rest.
    # Void-ratio ranges conditioned on soil type (anchors: Das Table 3.1 —
    # dense uniform sand e=0.45, loose e=0.8; stiff clay 0.6, soft 0.9-1.4).
    e_ranges = {
        "sand": (0.45, 0.85),
        "silt": (0.50, 0.90),
        "inorganic clay": (0.55, 0.95),
    }
    soil_type = random.choice(list(SPECIFIC_GRAVITY_RANGES.keys()))
    gs_lo, gs_hi = SPECIFIC_GRAVITY_RANGES[soil_type]
    Gs = round(random.uniform(gs_lo, gs_hi), 2)
    e_true = round(random.uniform(*e_ranges[soil_type]), 2)

    # Per-sample feasibility bounds on S so the presented w (= S*e/Gs) always
    # lands inside [5%, 33%] after rounding (small pre-rounding margins).
    s_lo = max(0.25, 0.051 * Gs / e_true)
    s_hi = min(0.92, 0.328 * Gs / e_true)
    S_true = random.uniform(s_lo, s_hi)

    w_true = S_true * e_true / Gs                # from S*e = w*Gs
    gamma_d_true = Gs * gamma_w / (1 + e_true)
    gamma_true = gamma_d_true * (1 + w_true)

    # Presented (rounded) given values — the question states exactly these.
    w_pct = round(w_true * 100, 1)               # moisture content, %
    gamma = round(gamma_true, 2)                 # moist unit weight, kN/m^3

    # 2. Core computation — round-then-recompute at EVERY step: each value in
    # the chain derives from the previously displayed (rounded) value, so the
    # printed arithmetic reproduces exactly.
    w = round(w_pct / 100.0, 3)
    gamma_d = round(gamma / (1 + w), 2)
    e = round((Gs * gamma_w / gamma_d) - 1, 3)
    S_frac = round(w * Gs / e, 3)
    S_pct = round(S_frac * 100, 1)

    # Physical bounds (docstring) enforced on the presented/recomputed chain.
    assert 2.60 <= Gs <= 2.80, f"Gs out of bounds: {Gs}"
    assert 14.0 <= gamma <= 23.0, f"moist unit weight out of bounds: {gamma}"
    assert 5.0 <= w_pct <= 33.0, f"moisture content out of bounds: {w_pct}"
    assert 0.30 <= e <= 1.10, f"recomputed void ratio out of bounds: {e}"
    assert 15.0 <= S_pct <= 100.0, f"recomputed saturation out of bounds: {S_pct}"

    # 3. Serialize question and gold trace.
    question = (
        f"A moist sample of {soil_type} taken from a borrow area has a total "
        f"(moist) unit weight of {gamma:.2f} kN/m^3 and a moisture content of "
        f"{w_pct:.1f}%. The specific gravity of the soil solids is {Gs:.2f}. "
        f"Taking the unit weight of water as {gamma_w:.2f} kN/m^3, determine "
        f"the degree of saturation of the sample in percent."
    )

    solution = (
        f"**Given:**\n"
        f"Moist unit weight (gamma): {gamma:.2f} kN/m^3\n"
        f"Moisture content (w): {w_pct:.1f}% = {w:.3f}\n"
        f"Specific gravity of solids (Gs): {Gs:.2f}\n"
        f"Unit weight of water (gamma_w): {gamma_w:.2f} kN/m^3\n\n"
        f"**Step 1:** Compute the dry unit weight.\n"
        f"The dry unit weight is obtained by dividing the total unit weight "
        f"by (1 + w):\n"
        f"gamma_d = gamma / (1 + w) = {gamma:.2f} / (1 + {w:.3f}) "
        f"= {gamma_d:.2f} kN/m^3\n\n"
        f"**Step 2:** Compute the void ratio.\n"
        f"From gamma_d = Gs * gamma_w / (1 + e), solving for e:\n"
        f"e = (Gs * gamma_w / gamma_d) - 1 "
        f"= ({Gs:.2f} * {gamma_w:.2f} / {gamma_d:.2f}) - 1 = {e:.3f}\n\n"
        f"**Step 3:** Compute the degree of saturation.\n"
        f"Using the relation S * e = w * Gs:\n"
        f"S = w * Gs / e = ({w:.3f} * {Gs:.2f}) / {e:.3f} "
        f"= {S_frac:.3f} = {S_pct:.1f}%\n\n"
        f"**Answer:** The degree of saturation is {S_pct:.1f} %"
    )

    return question, solution


# Template 2 (Easy) — Area B1: Phase Relationships & Index Properties
def template_relative_density_of_sand():
    """
    Relative Density of a Natural Sand Deposit

    Scenario:
        A clean sand deposit has a measured field dry unit weight, and
        laboratory tests provide the maximum and minimum void ratios. The
        field void ratio follows from the dry unit weight,

            e = (Gs * gamma_w / gamma_d) - 1

        and the relative density is

            Dr = (e_max - e) / (e_max - e_min)

    Difficulty: Easy
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Section 3.7 (Relative Density). Applies to cohesionless soils only,
        so the template samples sand exclusively; Gs per
        SPECIFIC_GRAVITY_RANGES["sand"].
    Physical bounds: e_min in [0.40, 0.52], e_max = e_min + [0.28, 0.42]
        (so e_max <= 0.94); target Dr sampled in [25%, 88%] so the
        recomputed field void ratio always lies strictly between e_min and
        e_max; presented dry unit weight in [13, 19] kN/m^3.

    Returns:
        tuple: (question, solution)
    """
    gamma_w = UNIT_WEIGHT_WATER_KN_M3

    # 1. Parameterize (sand only — relative density is defined for
    # cohesionless soils). e_min/e_max anchored to Das Table 3.1 (dense
    # uniform sand e = 0.45, loose e = 0.8).
    gs_lo, gs_hi = SPECIFIC_GRAVITY_RANGES["sand"]
    Gs = round(random.uniform(gs_lo, gs_hi), 2)
    e_min = round(random.uniform(0.40, 0.52), 2)
    e_max = round(e_min + random.uniform(0.28, 0.42), 2)
    Dr_true = random.uniform(0.25, 0.88)
    e_true = e_max - Dr_true * (e_max - e_min)

    gamma_d_true = Gs * gamma_w / (1 + e_true)
    gamma_d = round(gamma_d_true, 2)     # presented field dry unit weight

    # 2. Core computation — round-then-recompute at every step. e carries
    # 4 decimals: the small (e_max - e_min) denominator amplifies e-rounding
    # error by a factor of ~2.4-3.6, so 3 decimals would let the gold answer
    # drift beyond recomputation tolerance (R2 finding, cycle 1).
    e = round((Gs * gamma_w / gamma_d) - 1, 4)
    Dr_frac = round((e_max - e) / (e_max - e_min), 4)
    Dr_pct = round(Dr_frac * 100, 2)

    assert 0.40 <= e_min < e < e_max <= 0.94, (
        f"field void ratio outside (e_min, e_max): {e_min}, {e}, {e_max}")
    assert 13.0 <= gamma_d <= 19.0, f"dry unit weight out of bounds: {gamma_d}"
    assert 15.0 <= Dr_pct <= 92.0, f"relative density out of bounds: {Dr_pct}"

    # 3. Serialize.
    question = (
        f"A natural deposit of clean sand has a field dry unit weight of "
        f"{gamma_d:.2f} kN/m^3. Laboratory tests on the same sand give a "
        f"maximum void ratio of {e_max:.2f} and a minimum void ratio of "
        f"{e_min:.2f}. The specific gravity of the soil solids is {Gs:.2f} "
        f"and the unit weight of water is {gamma_w:.2f} kN/m^3. Determine "
        f"the relative density of the deposit in percent."
    )

    solution = (
        f"**Given:**\n"
        f"Field dry unit weight (gamma_d): {gamma_d:.2f} kN/m^3\n"
        f"Maximum void ratio (e_max): {e_max:.2f}\n"
        f"Minimum void ratio (e_min): {e_min:.2f}\n"
        f"Specific gravity of solids (Gs): {Gs:.2f}\n"
        f"Unit weight of water (gamma_w): {gamma_w:.2f} kN/m^3\n\n"
        f"**Step 1:** Compute the field void ratio from the dry unit weight.\n"
        f"From gamma_d = Gs * gamma_w / (1 + e), solving for e:\n"
        f"e = (Gs * gamma_w / gamma_d) - 1 "
        f"= ({Gs:.2f} * {gamma_w:.2f} / {gamma_d:.2f}) - 1 = {e:.4f}\n\n"
        f"**Step 2:** Compute the relative density.\n"
        f"Dr = (e_max - e) / (e_max - e_min) "
        f"= ({e_max:.2f} - {e:.4f}) / ({e_max:.2f} - {e_min:.2f}) "
        f"= {Dr_frac:.4f}\n\n"
        f"**Step 3:** Express the relative density as a percentage.\n"
        f"Dr = {Dr_frac:.4f} * 100 = {Dr_pct:.2f}%\n\n"
        f"**Answer:** The relative density is {Dr_pct:.2f} %"
    )

    return question, solution


# Template 3 (Intermediate) — Area B1: Phase Relationships & Index Properties
def template_borrow_pit_fill_volume():
    """
    Borrow-Pit Volume for a Compacted Fill

    Scenario:
        A compacted fill of specified volume and target void ratio must be
        constructed from borrow-area soil whose moist unit weight and
        moisture content are known. Because the weight of soil solids is
        conserved between the borrow pit and the fill, the chain is:

            gamma_d,fill   = Gs * gamma_w / (1 + e_fill)
            gamma_d,borrow = gamma_borrow / (1 + w_borrow)
            W_s            = gamma_d,fill * V_fill
            V_borrow       = W_s / gamma_d,borrow

    Difficulty: Intermediate
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Ch. 3 weight-volume relationships applied across two states
        (borrow-to-fill earthwork typology of Ch. 3/Ch. 6 problems).
    Physical bounds: e_fill in [0.40, 0.60] for sand, [0.45, 0.60] for silt;
        e_borrow = e_fill + [0.18, 0.45], capped at 0.90 (sand) / 1.05
        (silt) so the borrow deposit stays physically representable (borrow
        always looser than fill, so V_borrow > V_fill by construction);
        w_borrow in [8%, 18%]; presented unit weights in [12, 22] kN/m^3;
        1.05 <= V_borrow / V_fill <= 1.60.

    Returns:
        tuple: (question, solution)
    """
    gamma_w = UNIT_WEIGHT_WATER_KN_M3

    # 1. Parameterize: granular fill soils (sand or silt).
    soil_type = random.choice(["sand", "silt"])
    gs_lo, gs_hi = SPECIFIC_GRAVITY_RANGES[soil_type]
    Gs = round(random.uniform(gs_lo, gs_hi), 2)
    # Soil-conditioned ranges (R1, cycle 1): a natural sand deposit cannot be
    # looser than its e_max (~0.9), and a silt compacted to e = 0.40 would
    # demand near-modified-Proctor density — cap/floor per soil type.
    e_fill_ranges = {"sand": (0.40, 0.60), "silt": (0.45, 0.60)}
    e_borrow_caps = {"sand": 0.90, "silt": 1.05}
    e_fill = round(random.uniform(*e_fill_ranges[soil_type]), 2)
    # Sample the looseness spread up to the cap (no min() clamp: clamping
    # concentrated ~21% of sand draws exactly at the cap — R1, cycle 2).
    spread_hi = min(0.45, e_borrow_caps[soil_type] - e_fill)
    e_borrow = round(e_fill + random.uniform(0.18, spread_hi), 2)
    w_pct = round(random.uniform(8.0, 18.0), 1)
    V_fill = random.randint(60, 400) * 10          # fill volume, m^3

    # Input-form branch (Stage D remediation): the borrow soil is described
    # either by its moist state (gamma_borrow, w) or by its phase state
    # (e_borrow) — the sub-chain that reaches gamma_d,borrow changes with
    # the form, so the reasoning path is parameter-dependent.
    data_form = random.choice(["moist", "phase"])
    gamma_d_borrow_true = Gs * gamma_w / (1 + e_borrow)
    gamma_borrow = round(gamma_d_borrow_true * (1 + w_pct / 100.0), 2)

    # 2. Core computation — round-then-recompute at every step.
    w = round(w_pct / 100.0, 3)
    gamma_d_fill = round(Gs * gamma_w / (1 + e_fill), 2)
    if data_form == "moist":
        gamma_d_borrow = round(gamma_borrow / (1 + w), 2)
    else:
        gamma_d_borrow = round(Gs * gamma_w / (1 + e_borrow), 2)
    W_s = round(gamma_d_fill * V_fill, 1)
    V_borrow = round(W_s / gamma_d_borrow, 0)

    assert 12.0 <= gamma_d_borrow < gamma_d_fill <= 22.0, (
        f"unit-weight ordering violated: {gamma_d_borrow}, {gamma_d_fill}")
    assert 8.0 <= w_pct <= 18.0, f"moisture content out of bounds: {w_pct}"
    assert 1.05 <= V_borrow / V_fill <= 1.60, (
        f"borrow/fill volume ratio implausible: {V_borrow / V_fill}")

    # 3. Serialize.
    if data_form == "moist":
        borrow_text = (
            f"has a moist unit weight of {gamma_borrow:.2f} kN/m^3 and a "
            f"moisture content of {w_pct:.1f}%")
        given_borrow = (
            f"Borrow moist unit weight (gamma_borrow): {gamma_borrow:.2f} "
            f"kN/m^3\n"
            f"Borrow moisture content (w): {w_pct:.1f}% = {w:.3f}\n")
        step2 = (
            f"**Step 2:** Compute the dry unit weight of the borrow soil "
            f"from its moist state.\n"
            f"gamma_d,borrow = gamma_borrow / (1 + w) "
            f"= {gamma_borrow:.2f} / (1 + {w:.3f}) "
            f"= {gamma_d_borrow:.2f} kN/m^3\n\n")
    else:
        borrow_text = f"has an in-situ void ratio of {e_borrow:.2f}"
        given_borrow = f"Borrow void ratio (e_borrow): {e_borrow:.2f}\n"
        step2 = (
            f"**Step 2:** Compute the dry unit weight of the borrow soil "
            f"from its void ratio.\n"
            f"gamma_d,borrow = Gs * gamma_w / (1 + e_borrow) "
            f"= {Gs:.2f} * {gamma_w:.2f} / (1 + {e_borrow:.2f}) "
            f"= {gamma_d_borrow:.2f} kN/m^3\n\n")

    question = (
        f"A highway embankment requires {V_fill} m^3 of compacted fill with "
        f"a target void ratio of {e_fill:.2f}. The fill will be constructed "
        f"from a borrow area whose {soil_type} {borrow_text}. "
        f"The specific gravity of the soil solids is {Gs:.2f} and the unit "
        f"weight of water is {gamma_w:.2f} kN/m^3. Determine the volume of "
        f"borrow soil, in cubic meters, that must be excavated to build "
        f"the fill."
    )

    solution = (
        f"**Given:**\n"
        f"Required fill volume (V_fill): {V_fill} m^3\n"
        f"Target fill void ratio (e_fill): {e_fill:.2f}\n"
        f"{given_borrow}"
        f"Specific gravity of solids (Gs): {Gs:.2f}\n"
        f"Unit weight of water (gamma_w): {gamma_w:.2f} kN/m^3\n\n"
        f"**Step 1:** Compute the required dry unit weight of the compacted "
        f"fill.\n"
        f"gamma_d,fill = Gs * gamma_w / (1 + e_fill) "
        f"= {Gs:.2f} * {gamma_w:.2f} / (1 + {e_fill:.2f}) "
        f"= {gamma_d_fill:.2f} kN/m^3\n\n"
        f"{step2}"
        f"**Step 3:** Compute the weight of solids required in the fill.\n"
        f"The weight of soil solids is conserved from borrow pit to fill:\n"
        f"W_s = gamma_d,fill * V_fill = {gamma_d_fill:.2f} * {V_fill} "
        f"= {W_s:.1f} kN\n\n"
        f"**Step 4:** Compute the required borrow volume.\n"
        f"V_borrow = W_s / gamma_d,borrow = {W_s:.1f} / {gamma_d_borrow:.2f} "
        f"= {V_borrow:.0f} m^3\n\n"
        f"**Answer:** The required borrow volume is {V_borrow:.0f} m^3"
    )

    return question, solution
