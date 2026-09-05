import math
import random

from pilot.templates.branches.civil_engineering.constants import (
    SKEMPTON_CC_COEFF,
    SKEMPTON_CC_OFFSET,
    UNIT_WEIGHT_WATER_KN_M3,
)


# Template 7 (Advanced) — Area B3: Stress Distribution & Consolidation
def template_primary_consolidation_settlement():
    """
    Primary Consolidation Settlement Under a Raft Foundation (2:1 Method)

    Scenario:
        A normally consolidated clay layer lies beneath moist and saturated
        sand layers. A rectangular raft at the ground surface applies a
        uniform net pressure. The stress increase at the clay mid-depth is
        obtained by the 2:1 distribution method before the settlement
        computation:

            Cc          = 0.009 * (LL - 10)                    (Skempton)
            sigma'0     = gamma_m*z1 + gamma'_sand*z2 + gamma'_clay*Hc/2
            z           = z1 + z2 + Hc/2
            delta_sigma = q0*B*L / ((B + z)*(L + z))           (2:1 method)
            Sc          = Cc*Hc/(1+e0) * log10((sigma'0+delta_sigma)/sigma'0)

    Difficulty: Advanced
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Ch. 10 (stress increase, 2:1 approximate method), Ch. 11 (primary
        consolidation settlement of NC clay), Ch. 9 (effective stress);
        Skempton correlation per constants.SKEMPTON_CC_COEFF.
    Physical bounds: LL in [35, 65]; e0 coupled to LL so the implied
        saturated water content w = e0/Gs stays between 0.55*LL and
        0.95*LL (liquidity consistent with an NC clay, never above the
        liquid limit); raft B in [6, 12] m, L = B + [2, 6] m; q0 derived
        with per-sample feasibility so q0 in [80, 500] kPa, the stress
        ratio in [1.30, 2.00] (capped by the undrained-capacity screen:
        delta_sigma <= 1.0*sigma'0 < ~1.1*sigma'0 capacity), and the
        settlement in [30, 400] mm (asserts hold the same bands with
        rounding margins; the w/LL band edges may drift ~0.005 from 2-dp
        rounding of e0).

    Returns:
        tuple: (question, solution)
    """
    gamma_w = UNIT_WEIGHT_WATER_KN_M3

    # 1. Parameterize the profile.
    z1 = round(random.uniform(1.5, 3.0), 1)        # moist sand above WT, m
    z2 = round(random.uniform(2.0, 5.0), 1)        # saturated sand, m
    Hc = round(random.uniform(2.0, 5.0), 1)        # clay thickness, m
    gamma_m = round(random.uniform(16.5, 18.0), 1)     # moist sand, kN/m^3
    gamma_sat_sand = round(random.uniform(19.0, 20.5), 1)
    LL = random.randint(35, 65)                    # liquid limit, %
    Gs_clay = round(random.uniform(2.70, 2.80), 2)
    # Couple e0 to LL (R1, cycle 1): the saturated water content implied by
    # e0 (w = e0/Gs) must sit inside [0.55*LL, 0.95*LL] so the clay is a
    # plausible NC deposit (liquidity index < 1, not quick; not
    # atypically stiff either). Window is non-empty across LL in [35, 65]
    # and Gs in [2.70, 2.80].
    e0_lo = max(0.80, 0.55 * (LL / 100.0) * Gs_clay)
    e0_hi = min(1.20, 0.95 * (LL / 100.0) * Gs_clay)
    e0 = round(random.uniform(e0_lo, e0_hi), 2)
    # Clay saturated unit weight derived from its own phase state.
    gamma_sat_clay = round((Gs_clay + e0) * gamma_w / (1 + e0), 2)

    # Raft geometry.
    B = round(random.uniform(6.0, 12.0), 1)
    L = round(B + random.uniform(2.0, 6.0), 1)

    # Initial effective stress at clay mid-depth, from presented values.
    gamma_sub_sand = round(gamma_sat_sand - gamma_w, 2)
    gamma_sub_clay = round(gamma_sat_clay - gamma_w, 2)
    sigma0 = round(
        gamma_m * z1 + gamma_sub_sand * z2 + gamma_sub_clay * Hc / 2, 1)
    z = round(z1 + z2 + Hc / 2, 2)
    F = B * L / ((B + z) * (L + z))                # 2:1 attenuation factor

    # Per-sample feasibility for the raft pressure (lessons 1, 9): sample
    # log10(stress ratio) inside the window that jointly keeps the ratio in
    # [1.30, 3.50], the settlement in [0.03, 0.40] m, and q0 in
    # [80, 500] kPa.
    Cc = round(SKEMPTON_CC_COEFF * (LL - SKEMPTON_CC_OFFSET), 3)
    M = Cc * Hc / (1 + e0)
    # Ratio capped at 2.00 (R1, cycle 2): delta_sigma <= 1.0*sigma'0 keeps
    # the increment below the NC clay's undrained capacity screen
    # (~(1.1-1.6)*sigma'0 with su ~= 0.22-0.25*sigma'v0), so the clay
    # consolidates rather than failing in shear.
    log_r_lo = max(math.log10(1.30), 0.03 / M,
                   math.log10(1 + 80.0 * F / sigma0))
    log_r_hi = min(math.log10(2.00), 0.40 / M,
                   math.log10(1 + 500.0 * F / sigma0))
    assert log_r_lo < log_r_hi, (
        f"empty feasibility window: {log_r_lo}, {log_r_hi}")
    log_r = random.uniform(log_r_lo, log_r_hi)
    q0 = round(sigma0 * (10 ** log_r - 1) / F, 0)  # presented raft pressure

    # 2. Core computation — round-then-recompute at every step.
    delta_sigma = round(q0 * B * L / ((B + z) * (L + z)), 1)
    sigma_f = round(sigma0 + delta_sigma, 1)
    ratio = sigma_f / sigma0
    log_term = round(math.log10(sigma_f / sigma0), 4)
    Sc_m = round(Cc * Hc / (1 + e0) * log_term, 4)
    Sc_mm = round(Sc_m * 1000, 0)

    assert 30.0 <= sigma0 <= 150.0, f"sigma'0 out of bounds: {sigma0}"
    assert 75.0 <= q0 <= 505.0, f"raft pressure out of bounds: {q0}"
    assert 1.28 <= ratio <= 2.05, f"stress ratio out of bounds: {ratio}"
    assert 28.0 <= Sc_mm <= 405.0, f"settlement out of bounds: {Sc_mm}"

    # 3. Serialize.
    question = (
        f"A soil profile consists of, from the ground surface down: "
        f"{z1:.1f} m of moist sand (unit weight {gamma_m:.1f} kN/m^3) above "
        f"the water table; {z2:.1f} m of saturated sand (saturated unit "
        f"weight {gamma_sat_sand:.1f} kN/m^3); and a {Hc:.1f} m thick layer "
        f"of normally consolidated clay (saturated unit weight "
        f"{gamma_sat_clay:.2f} kN/m^3, initial void ratio {e0:.2f}, liquid "
        f"limit {LL}%). The water table is at the bottom of the moist sand "
        f"layer. A rectangular raft foundation, {B:.1f} m by {L:.1f} m in "
        f"plan, applies a uniform net pressure of {q0:.0f} kPa at the "
        f"ground surface. Using the 2:1 stress distribution method for the "
        f"stress increase at the middle of the clay layer, the Skempton "
        f"correlation Cc = 0.009(LL - 10), the settlement relation "
        f"Sc = [Cc*Hc/(1+e0)] * log10(sigma'f/sigma'0), and a unit weight "
        f"of water of {gamma_w:.2f} kN/m^3, determine the primary "
        f"consolidation settlement of the clay layer in millimeters."
    )

    solution = (
        f"**Given:**\n"
        f"Moist sand: z1 = {z1:.1f} m, gamma = {gamma_m:.1f} kN/m^3\n"
        f"Saturated sand: z2 = {z2:.1f} m, gamma_sat = "
        f"{gamma_sat_sand:.1f} kN/m^3\n"
        f"Clay: Hc = {Hc:.1f} m, gamma_sat = {gamma_sat_clay:.2f} kN/m^3, "
        f"e0 = {e0:.2f}, LL = {LL}%\n"
        f"Raft: {B:.1f} m x {L:.1f} m, net pressure q0 = {q0:.0f} kPa at "
        f"the surface\n"
        f"Unit weight of water (gamma_w): {gamma_w:.2f} kN/m^3\n\n"
        f"**Step 1:** Estimate the compression index from the Skempton "
        f"correlation.\n"
        f"Cc = 0.009 * (LL - 10) = 0.009 * ({LL} - 10) = {Cc:.3f}\n\n"
        f"**Step 2:** Compute the submerged unit weights below the water "
        f"table.\n"
        f"Sand: gamma'_sand = {gamma_sat_sand:.1f} - {gamma_w:.2f} = "
        f"{gamma_sub_sand:.2f} kN/m^3\n"
        f"Clay: gamma'_clay = {gamma_sat_clay:.2f} - {gamma_w:.2f} = "
        f"{gamma_sub_clay:.2f} kN/m^3\n\n"
        f"**Step 3:** Compute the initial effective stress at the middle of "
        f"the clay layer.\n"
        f"sigma'0 = {gamma_m:.1f} * {z1:.1f} + {gamma_sub_sand:.2f} * "
        f"{z2:.1f} + {gamma_sub_clay:.2f} * {Hc / 2:.2f} = {sigma0:.1f} "
        f"kPa\n\n"
        f"**Step 4:** Compute the depth from the raft to the middle of the "
        f"clay layer.\n"
        f"z = z1 + z2 + Hc/2 = {z1:.1f} + {z2:.1f} + {Hc / 2:.2f} = "
        f"{z:.2f} m\n\n"
        f"**Step 5:** Compute the stress increase at that depth by the 2:1 "
        f"method.\n"
        f"The load spreads over an area (B + z) by (L + z):\n"
        f"delta_sigma = q0 * B * L / ((B + z) * (L + z)) "
        f"= {q0:.0f} * {B:.1f} * {L:.1f} / (({B:.1f} + {z:.2f}) * "
        f"({L:.1f} + {z:.2f})) = {delta_sigma:.1f} kPa\n\n"
        f"**Step 6:** Compute the final effective stress at the clay "
        f"mid-depth.\n"
        f"sigma'f = sigma'0 + delta_sigma = {sigma0:.1f} + "
        f"{delta_sigma:.1f} = {sigma_f:.1f} kPa\n\n"
        f"**Step 7:** Compute the logarithmic stress-ratio term.\n"
        f"log10(sigma'f / sigma'0) = log10({sigma_f:.1f} / {sigma0:.1f}) = "
        f"{log_term:.4f}\n\n"
        f"**Step 8:** Compute the primary consolidation settlement.\n"
        f"Sc = Cc * Hc / (1 + e0) * log10(sigma'f/sigma'0) "
        f"= {Cc:.3f} * {Hc:.1f} / (1 + {e0:.2f}) * {log_term:.4f} "
        f"= {Sc_m:.4f} m\n\n"
        f"**Step 9:** Convert to millimeters.\n"
        f"Sc = {Sc_m:.4f} * 1000 = {Sc_mm:.0f} mm\n\n"
        f"**Answer:** The primary consolidation settlement is "
        f"{Sc_mm:.0f} mm"
    )

    return question, solution


# Template 8 (Advanced) — Area B3: Stress Distribution & Consolidation
def template_time_rate_of_consolidation():
    """
    Field Consolidation Time from a Laboratory Oedometer Test

    Scenario:
        A doubly drained laboratory specimen of the field clay reaches 50%
        consolidation in a measured time. Because the time factor
        Tv = cv*t / H_dr^2 is the same dimensionless group in the lab and
        the field, the field time for a target degree of consolidation
        follows by similitude without computing cv explicitly:

            t_field = t50_lab * (Tv_field / T50) * (H_dr,field / H_dr,lab)^2

        with the time factors from the standard relations (stated in the
        question):
            U <= 60%:  Tv = (pi/4) * (U/100)^2
            U  > 60%:  Tv = 1.781 - 0.933 * log10(100 - U)

    Difficulty: Advanced
    Grounding: Das & Sobhan, Principles of Geotechnical Engineering, 9th ed.,
        Ch. 11 (time rate of consolidation; time-factor relations;
        lab-to-field extrapolation typology).
    Physical bounds: lab specimen 25 mm thick, doubly drained
        (H_dr,lab = 12.5 mm); field clay Hc in [2.0, 6.0] m with sampled
        drainage condition (double -> Hc/2, single -> Hc); U in
        {40, 55, 65, 80, 90}% (two of five values exercise the parabolic
        relation; 50 excluded as degenerate against the lab reference
        degree); lab t50 sampled with per-sample feasibility inside
        [8, 40] min so the field time stays in [0.08, 22] years (window
        provably non-empty: the similitude factor R lies in [4096, 995328]
        and the t50 window closes only outside [1072, 1.42e6]).

    Returns:
        tuple: (question, solution)
    """
    # 1. Parameterize; drainage condition and U each change the reasoning
    # path (drainage path AND which Tv relation applies).
    Hc = round(random.uniform(2.0, 6.0), 1)
    drainage = random.choice(["double", "single"])
    # 50 is excluded: it equals the lab reference degree, making the
    # time-factor ratio exactly 1 and bypassing the Tv-relation reasoning
    # (R2+R3, cycle 2).
    U = random.choice([40, 55, 65, 80, 90])
    H_lab_mm = 12.5                                 # 25 mm specimen, doubly drained
    T50 = round(math.pi / 4.0 * 0.25, 4)            # = 0.1963

    H_dr_true = Hc / 2.0 if drainage == "double" else Hc
    if U <= 60:
        Tv_true = math.pi / 4.0 * (U / 100.0) ** 2
    else:
        Tv_true = 1.781 - 0.933 * math.log10(100 - U)

    # Per-sample feasibility for the lab t50 (lesson 1): keep the field
    # time inside [0.08, 22] years given the similitude factor.
    ratio_sq_true = (H_dr_true * 1000.0 / H_lab_mm) ** 2
    R = (Tv_true / T50) * ratio_sq_true
    minutes_per_year = 525600.0                     # 60 * 24 * 365
    t50_lo = max(8.0, 1.02 * 0.08 * minutes_per_year / R)
    t50_hi = min(40.0, 0.98 * 22.0 * minutes_per_year / R)
    t50 = round(random.uniform(t50_lo, t50_hi), 0)

    # 2. Core computation — round-then-recompute at every step.
    H_dr = round(H_dr_true, 2)                      # field drainage path, m
    Tv = round(Tv_true, 4)
    ratio = round(H_dr * 1000.0 / H_lab_mm, 1)      # drainage-path ratio
    t_min = round(t50 * (Tv / T50) * ratio ** 2, 0)
    # 3 decimals below one year (Stage D audit: 2-dp quantization reached
    # ~1-6% relative on the shortest times).
    t_prec = 3 if t_min / minutes_per_year < 1.0 else 2
    t_years = round(t_min / minutes_per_year, t_prec)

    assert 0.12 <= Tv <= 0.90, f"time factor out of bounds: {Tv}"
    assert 0.05 <= t_years <= 23.0, f"field time out of bounds: {t_years}"

    drainage_text = (
        "sand layers above and below the clay (double drainage)"
        if drainage == "double"
        else "a sand layer above the clay and impervious rock below it "
             "(single drainage)"
    )

    if drainage == "double":
        step2 = (
            f"**Step 2:** Determine the field drainage path.\n"
            f"With permeable sand above and below, the clay drains from "
            f"both faces:\n"
            f"H_dr,field = Hc / 2 = {Hc:.1f} / 2 = {H_dr:.2f} m\n\n"
        )
    else:
        step2 = (
            f"**Step 2:** Determine the field drainage path.\n"
            f"With impervious rock below, the clay drains through its top "
            f"face only:\n"
            f"H_dr,field = Hc = {H_dr:.2f} m\n\n"
        )
    if U <= 60:
        step3 = (
            f"**Step 3:** Compute the field time factor for U = {U}%.\n"
            f"Since U <= 60%, Tv = (pi/4) * (U/100)^2 "
            f"= (pi/4) * ({U / 100.0:.2f})^2 = {Tv:.4f}\n\n"
        )
    else:
        step3 = (
            f"**Step 3:** Compute the field time factor for U = {U}%.\n"
            f"Since U > 60%, Tv = 1.781 - 0.933 * log10(100 - U) "
            f"= 1.781 - 0.933 * log10({100 - U}) = {Tv:.4f}\n\n"
        )

    question = (
        f"A {Hc:.1f} m thick saturated clay layer in the field is bounded "
        f"by {drainage_text}. In the laboratory, a 25 mm thick specimen of "
        f"the same clay, drained on both faces, reaches 50% average "
        f"consolidation in {t50:.0f} minutes. The time factor is defined "
        f"as Tv = cv*t / H_dr^2, where the specimen and the field layer "
        f"share the same coefficient of consolidation cv. Using the "
        f"time-factor relations Tv = (pi/4)(U/100)^2 for U <= 60% and "
        f"Tv = 1.781 - 0.933*log10(100 - U) for U > 60%, determine the "
        f"time in years required for the field layer to reach an average "
        f"degree of consolidation of {U}%. (Use 525,600 minutes per year.)"
    )

    solution = (
        f"**Given:**\n"
        f"Field clay thickness (Hc): {Hc:.1f} m, {drainage} drainage\n"
        f"Lab specimen: 25 mm thick, drained on both faces, t50 = "
        f"{t50:.0f} min\n"
        f"Target field degree of consolidation (U): {U}%\n\n"
        f"**Step 1:** Determine the laboratory drainage path and its time "
        f"factor.\n"
        f"The doubly drained 25 mm specimen has H_dr,lab = 25 / 2 = "
        f"{H_lab_mm:.1f} mm. At 50% consolidation (U <= 60%):\n"
        f"T50 = (pi/4) * (0.50)^2 = {T50:.4f}\n\n"
        f"{step2}"
        f"{step3}"
        f"**Step 4:** Form the drainage-path ratio between field and lab.\n"
        f"H_dr,field / H_dr,lab = {H_dr:.2f} m / {H_lab_mm:.1f} mm "
        f"= {H_dr * 1000:.0f} mm / {H_lab_mm:.1f} mm = {ratio:.1f}\n\n"
        f"**Step 5:** Apply similitude of the time factor "
        f"(Tv = cv*t/H_dr^2, same cv in lab and field).\n"
        f"t_field = t50 * (Tv / T50) * (H_dr,field / H_dr,lab)^2 "
        f"= {t50:.0f} * ({Tv:.4f} / {T50:.4f}) * ({ratio:.1f})^2 "
        f"= {t_min:.0f} min\n\n"
        f"**Step 6:** Convert to years.\n"
        f"t_field = {t_min:.0f} / 525600 = {t_years:.{t_prec}f} years\n\n"
        f"**Answer:** The required field consolidation time is "
        f"{t_years:.{t_prec}f} years"
    )

    return question, solution
