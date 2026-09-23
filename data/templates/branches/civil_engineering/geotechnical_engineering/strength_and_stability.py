import math
import random
from decimal import Decimal, ROUND_HALF_UP

from data.templates.branches.civil_engineering.constants import (
    FRICTION_ANGLE_RANGES_DEG,
    TERZAGHI_BEARING_FACTORS,
    TERZAGHI_MODIFIED_FACTORS,
)


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


# Template 9 (Easy) — Area B4: Shear Strength & Stability Applications
def template_infinite_slope_factor_of_safety():
    """
    Factor of Safety of a Dry Cohesionless Infinite Slope

    Scenario:
        A long, planar slope in dry cohesionless soil is analyzed with the
        infinite-slope model. With no cohesion and no seepage, the factor
        of safety reduces to

            FS = tan(phi') / tan(beta)

        independent of the failure-plane depth and unit weight.

    Difficulty: Easy
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Ch. 15 (slope stability; infinite slope without seepage). Friction
        angles per constants.FRICTION_ANGLE_RANGES_DEG (Das-cited).
    Physical bounds: phi' sampled per soil type from the Das ranges
        (27-45 deg); slope angle derived from a target FS in [1.15, 2.20]
        so beta always lies in [13, 41] deg (joint sampling range) and
        strictly below phi' (stable slope, FS > 1).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize: granular soils only (cohesionless model).
    soil_options = [
        "sand, rounded, loose",
        "sand, rounded, dense",
        "sand, angular, loose",
        "sand, angular, dense",
    ]
    soil_type = random.choice(soil_options)
    phi_lo, phi_hi = FRICTION_ANGLE_RANGES_DEG[soil_type]
    phi = random.randint(phi_lo, phi_hi)
    FS_true = random.uniform(1.15, 2.20)
    beta = round(math.degrees(math.atan(math.tan(math.radians(phi)) / FS_true)), 1)

    # 2. Core computation — round-then-recompute from the presented angles.
    tan_phi = round(math.tan(math.radians(phi)), 4)
    tan_beta = round(math.tan(math.radians(beta)), 4)
    FS = round(tan_phi / tan_beta, 3)

    assert 10.0 <= beta < phi, f"slope angle implausible: {beta} vs phi {phi}"
    assert 1.05 <= FS <= 2.35, f"factor of safety out of bounds: {FS}"

    # 3. Serialize.
    base_soil = soil_type.split(",")[0].strip()               # "sand"
    descriptor = soil_type.split(",", 1)[1].strip()           # "rounded, dense"
    question = (
        f"A natural slope in dry, cohesionless {base_soil} "
        f"({descriptor}) is inclined at {beta:.1f} degrees to the "
        f"horizontal. The drained friction angle of the soil is {phi} "
        f"degrees, and there is no seepage. Using the infinite-slope "
        f"model, determine the factor of safety against sliding."
    )

    solution = (
        f"**Given:**\n"
        f"Slope angle (beta): {beta:.1f} deg\n"
        f"Drained friction angle (phi'): {phi} deg\n"
        f"Soil: dry, cohesionless (c' = 0), no seepage\n\n"
        f"**Step 1:** Select the governing relation.\n"
        f"For an infinite slope in dry cohesionless soil, the sliding "
        f"resistance is purely frictional, and the factor of safety is "
        f"independent of depth and unit weight:\n"
        f"FS = tan(phi') / tan(beta)\n\n"
        f"**Step 2:** Evaluate the two tangents.\n"
        f"tan(phi') = tan({phi} deg) = {tan_phi:.4f}\n"
        f"tan(beta) = tan({beta:.1f} deg) = {tan_beta:.4f}\n\n"
        f"**Step 3:** Compute the factor of safety.\n"
        f"FS = {tan_phi:.4f} / {tan_beta:.4f} = {FS:.3f}\n\n"
        f"**Answer:** The factor of safety against sliding is {FS:.3f}"
    )

    return question, solution


# Template 10 (Intermediate) — Area B4: Shear Strength & Stability Applications
def template_terzaghi_strip_footing_bearing():
    """
    Allowable Bearing Capacity of a Strip Footing (General vs Local Shear)

    Scenario:
        A strip footing rests on a soil whose density/compaction state
        determines the bearing failure mode. The solver must first decide
        the governing regime from the stated site data, then apply the
        matching Terzaghi equation:

            general shear: qu = c'*Nc + q*Nq + 0.5*gamma*B*Ngamma
            local shear:   qu = (2/3)c'*N'c + q*N'q + 0.5*gamma*B*N'gamma

        with q = gamma*Df and FS = 3 for the allowable pressure. The
        failure-mode decision changes the equation, the factor set, and
        (for c' > 0) the mobilized cohesion — parameter-dependent
        reasoning-path branching.

    Difficulty: Intermediate
    Grounding: Das & Sobhan PGE 9th ed. Ch. 16 Table 16.1 (general-shear
        factors, Kumbhojkar Ngamma) and Das & Sivakugan PFE 9th ed.
        Table 3.2 + Eq. (3.9) (local-shear modified factors and equation),
        both transcribed on-disk into constants; factor values and both
        equations stated in the question per the self-containment policy.
    Physical bounds: two scenario families x two modes — sand (c' = 0):
        general = dense (Dr 72-90%, phi' in {30, 35}, gamma 18.0-19.5),
        local = loose-to-medium (Dr 30-60%, phi' = 30, gamma 16.5-17.5);
        clayey fill (c' in {10, 15, 20} kPa, phi' in {20, 25}): general =
        well-compacted (96-99% standard Proctor, gamma 18.0-19.0), local =
        lightly compacted (88-93% Proctor, gamma 16.5-17.5). Decision
        thresholds stated in the question (Dr > 70% / >= 95% Proctor).
        B in [1.0, 2.5] m, Df in [0.5, 1.5] m; qu floors 250 kPa (general)
        / 95 kPa (local), ceiling 2600 kPa.

    Trace integrity (Layer 0, 2026-09-23):
        The general-shear cohesion term c'*Nc (an integer times a 2-dp
        table factor) is exact at 2 dp and was printed at 1 dp, where it
        sat on a half-way tie on every c' = 15 draw (15 * 25.13 = 376.95);
        it is now bound half-up at 2 dp and printed at 2 dp in Step 4 and
        in the Step 5 sum, and qu, the sum of that 2-dp term and two 1-dp
        terms, is bound and printed at 2 dp with it. The local-shear chain
        stays at 1 dp ((2/3)c'*N'c takes six values, none on a tie). The
        surcharge and width terms are 4- and 5-dp products printed at 1 dp,
        so a draw that lands either on a 1-dp half-way tie is resampled
        rather than rounded either way (D-016/D-037). q_all = qu/3 is a
        quotient exact at no fixed display, and at 1 dp it tied on a third
        of c' = 15 general-shear draws; on the owner's decision of
        2026-09-23 (D-044) the answer display was lengthened to 2 dp
        rather than those draws resampled. A 2-dp q_all cannot tie: with
        qu = n/100, a 2-dp half-way point needs n/3 = k + 1/2, i.e.
        2n = 3(2k + 1), which no integer n satisfies.

    Returns:
        tuple: (question, solution)
    """
    # Bounded redraw: a draw whose cohesion, surcharge or width term lands
    # on a half-way tie at its display has no defensible gold answer and is
    # rejected (D-016).
    for _attempt in range(200):
        scenario = random.choice(["sand", "clayey"])
        mode = random.choice(["general", "local"])

        if scenario == "sand":
            c = 0
            if mode == "general":
                soil_desc = "a deep deposit of dense sand"
                Dr = random.randint(72, 90)
                # phi' = 35 only: Das dense rounded sand is 35-38 deg; pairing
                # Dr 72-90% with phi' = 30 contradicted the constants table
                # (R1, branching cycle).
                phi = 35
                gamma = round(random.uniform(18.0, 19.5), 1)
            else:
                soil_desc = "a deep deposit of loose to medium sand"
                Dr = random.randint(30, 60)
                # phi' = 30 = top of Das loose rounded range (27-30), i.e. the
                # loose/medium boundary — consistent with the descriptor.
                phi = 30
                gamma = round(random.uniform(16.5, 17.5), 1)
            datum = f"Field tests indicate a relative density of {Dr}%."
            rule = ("Assume general shear failure governs when the relative "
                    "density exceeds 70%; otherwise local shear failure "
                    "governs.")
            decision_evidence = (
                f"the relative density is {Dr}%, which is "
                f"{'greater' if mode == 'general' else 'not greater'} than 70%")
        else:
            c = random.choice([10, 15, 20])
            phi = random.choice([20, 25])
            if mode == "general":
                soil_desc = "a well-compacted clayey fill"
                proctor = random.randint(96, 99)
                gamma = round(random.uniform(18.0, 19.0), 1)
            else:
                soil_desc = "a lightly compacted clayey fill"
                proctor = random.randint(88, 93)
                gamma = round(random.uniform(16.5, 17.5), 1)
            datum = (f"The fill was placed at {proctor}% of its standard "
                     f"Proctor maximum dry density.")
            rule = ("Assume general shear failure governs when the fill is "
                    "compacted to at least 95% of standard Proctor density; "
                    "otherwise local shear failure governs.")
            decision_evidence = (
                f"the fill is at {proctor}% of standard Proctor density, "
                f"{'at or above' if mode == 'general' else 'below'} the 95% "
                f"threshold")

        Nc, Nq, Ngamma = TERZAGHI_BEARING_FACTORS[phi]
        Ncp, Nqp, Ngp = TERZAGHI_MODIFIED_FACTORS[phi]
        B = round(random.uniform(1.0, 2.5), 1)
        Df = round(random.uniform(0.5, 1.5), 1)
        FS = 3

        # Core computation — round-then-recompute at every step, along the
        # branch selected by the failure-mode decision.
        q_sur = round(gamma * Df, 2)
        if mode == "general":
            c_used = c
            # An integer c' times a 2-dp table factor is exact at 2 dp; at
            # 1 dp it sat on a half-way tie on every c' = 15 draw
            # (15 * 25.13 = 376.95), so it is carried and printed at 2 dp.
            tc_dp = 2
            term_c_exact = c_used * Nc
            term_q_exact = q_sur * Nq
            term_g_exact = 0.5 * gamma * B * Ngamma
        else:
            c_used = round(2.0 / 3.0 * c, 2)
            # (2/3)c' * N'c takes one of six values, none on a 1-dp tie.
            tc_dp = 1
            term_c_exact = c_used * Ncp
            term_q_exact = q_sur * Nqp
            term_g_exact = 0.5 * gamma * B * Ngp
        # The surcharge and width terms are 4- and 5-dp products printed at
        # 1 dp; the rare draw that lands on a half-way tie is redrawn.
        if (_is_display_tie(term_c_exact, tc_dp)
                or _is_display_tie(term_q_exact, 1)
                or _is_display_tie(term_g_exact, 1)):
            continue                    # no defensible gold answer; redraw
        term_c = _hu(term_c_exact, tc_dp)
        term_q = round(term_q_exact, 1)
        term_g = round(term_g_exact, 1)
        # The sum of the three displayed terms is exact at tc_dp.
        qu = _hu(term_c + term_q + term_g, tc_dp)
        # qu / 3 is exact at no fixed display, but with qu = n/100 a 2-dp
        # half-way point needs 2n = 3(2k + 1), impossible for integer n; the
        # answer is printed at 2 dp, where it cannot tie (D-044).
        q_all = round(qu / FS, 2)
        break
    else:
        raise RuntimeError(
            "terzaghi_strip_footing_bearing: no closing sample in 200 draws")

    if mode == "general":
        assert 250.0 <= qu <= 2600.0, f"general-shear qu out of bounds: {qu}"
    else:
        assert 95.0 <= qu <= 900.0, f"local-shear qu out of bounds: {qu}"
    assert 30.0 <= q_all <= 900.0, f"allowable capacity out of bounds: {q_all}"

    question = (
        f"A strip footing of width B = {B:.1f} m is founded at a depth "
        f"Df = {Df:.1f} m in {soil_desc} with unit weight gamma = "
        f"{gamma:.1f} kN/m^3, effective cohesion c' = {c} kPa, and drained "
        f"friction angle phi' = {phi} degrees. {datum} The groundwater "
        f"table is well below the footing. Terzaghi's strip-footing "
        f"equations are "
        f"qu = c'*Nc + q*Nq + 0.5*gamma*B*Ngamma for general shear "
        f"failure and qu = (2/3)c'*N'c + q*N'q + 0.5*gamma*B*N'gamma for "
        f"local shear failure, where q = gamma*Df. For phi' = {phi} "
        f"degrees: Nc = {Nc:.2f}, Nq = {Nq:.2f}, Ngamma = {Ngamma:.2f}; "
        f"N'c = {Ncp:.2f}, N'q = {Nqp:.2f}, N'gamma = {Ngp:.2f}. {rule} "
        f"The allowable bearing capacity is q_all = qu / FS. Determine the "
        f"allowable bearing capacity in kPa for a factor of safety of "
        f"{FS}."
    )

    if mode == "general":
        step2 = (
            f"**Step 2:** Select the governing equation and factors.\n"
            f"qu = c'*Nc + q*Nq + 0.5*gamma*B*Ngamma, with Nc = {Nc:.2f}, "
            f"Nq = {Nq:.2f}, Ngamma = {Ngamma:.2f}. The full effective "
            f"cohesion c' = {c} kPa is mobilized.\n\n"
        )
        terms = (
            f"Cohesion term: c' * Nc = {c} * {Nc:.2f} = {term_c:.{tc_dp}f} kPa\n"
            f"Surcharge term: q * Nq = {q_sur:.2f} * {Nq:.2f} = "
            f"{term_q:.1f} kPa\n"
            f"Width term: 0.5 * gamma * B * Ngamma = 0.5 * {gamma:.1f} * "
            f"{B:.1f} * {Ngamma:.2f} = {term_g:.1f} kPa\n\n"
        )
    else:
        step2 = (
            f"**Step 2:** Select the governing equation and factors.\n"
            f"qu = (2/3)c'*N'c + q*N'q + 0.5*gamma*B*N'gamma, with "
            f"N'c = {Ncp:.2f}, N'q = {Nqp:.2f}, N'gamma = {Ngp:.2f}. The "
            f"mobilized cohesion is reduced: (2/3)c' = (2/3) * {c} = "
            f"{c_used:.2f} kPa.\n\n"
        )
        terms = (
            f"Cohesion term: (2/3)c' * N'c = {c_used:.2f} * {Ncp:.2f} = "
            f"{term_c:.{tc_dp}f} kPa\n"
            f"Surcharge term: q * N'q = {q_sur:.2f} * {Nqp:.2f} = "
            f"{term_q:.1f} kPa\n"
            f"Width term: 0.5 * gamma * B * N'gamma = 0.5 * {gamma:.1f} * "
            f"{B:.1f} * {Ngp:.2f} = {term_g:.1f} kPa\n\n"
        )

    solution = (
        f"**Given:**\n"
        f"Footing width (B): {B:.1f} m, depth (Df): {Df:.1f} m\n"
        f"Soil: gamma = {gamma:.1f} kN/m^3, c' = {c} kPa, phi' = {phi} deg\n"
        f"Site data: {datum}\n"
        f"General-shear factors: Nc = {Nc:.2f}, Nq = {Nq:.2f}, "
        f"Ngamma = {Ngamma:.2f}\n"
        f"Local-shear factors: N'c = {Ncp:.2f}, N'q = {Nqp:.2f}, "
        f"N'gamma = {Ngp:.2f}\n"
        f"Factor of safety (FS): {FS}\n\n"
        f"**Step 1:** Decide the failure mode.\n"
        f"Since {decision_evidence}, {mode} shear failure governs.\n\n"
        f"{step2}"
        f"**Step 3:** Compute the effective surcharge at the footing "
        f"base.\n"
        f"q = gamma * Df = {gamma:.1f} * {Df:.1f} = {q_sur:.2f} kPa\n\n"
        f"**Step 4:** Evaluate the three terms of the selected equation.\n"
        f"{terms}"
        f"**Step 5:** Sum the terms to obtain the ultimate bearing "
        f"capacity.\n"
        f"qu = {term_c:.{tc_dp}f} + {term_q:.1f} + {term_g:.1f} = "
        f"{qu:.{tc_dp}f} kPa\n\n"
        f"**Step 6:** Apply the factor of safety.\n"
        f"q_all = qu / FS = {qu:.{tc_dp}f} / {FS} = {q_all:.2f} kPa\n\n"
        f"**Answer:** The allowable bearing capacity is {q_all:.2f} kPa"
    )

    return question, solution
