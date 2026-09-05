import math
import random

from pilot.templates.branches.civil_engineering.constants import (
    RATIONAL_C,
    SCS_CURVE_NUMBERS,
)

# Land uses whose HEC-22 C ranges suit small urban catchments, with
# natural-language display phrases (R3, cycle 1: raw taxonomy labels with
# internal colons read awkwardly in prose).
_RATIONAL_USES = {
    "business: downtown": "downtown business district",
    "business: neighborhood": "neighborhood business district",
    "residential: single-family": "single-family residential area",
    "residential: multi-unit detached": "multi-unit detached residential area",
    "residential: apartments": "apartment residential area",
    "industrial: light": "light industrial area",
    "parks and cemeteries": "park and cemetery area",
    "playgrounds": "playground area",
}

# (description, HSG, CN) combos restricted to CN >= 61 so the sampled
# rainfall window is always comfortably above the initial abstraction.
_SCS_COMBOS = [
    ("commercial and business (85% impervious)", "B", 92),
    ("commercial and business (85% impervious)", "C", 94),
    ("residential: 1/4 acre lots (38% impervious)", "A", 61),
    ("residential: 1/4 acre lots (38% impervious)", "B", 75),
    ("residential: 1/4 acre lots (38% impervious)", "C", 83),
    ("residential: 1/8 acre lots (65% impervious)", "B", 85),
    ("open space, good condition (>75% grass)", "B", 61),
    ("open space, good condition (>75% grass)", "C", 74),
    ("industrial (72% impervious)", "B", 88),
    ("streets: gravel", "B", 85),
]

_SCS_DISPLAY = {
    "commercial and business (85% impervious)": "a commercial and business district (85% impervious)",
    "residential: 1/4 acre lots (38% impervious)": "a residential district of 1/4-acre lots (38% impervious)",
    "residential: 1/8 acre lots (65% impervious)": "a residential district of 1/8-acre lots (65% impervious)",
    "open space, good condition (>75% grass)": "open space in good condition (>75% grass cover)",
    "industrial (72% impervious)": "an industrial district (72% impervious)",
    "streets: gravel": "a gravel-road corridor",
}


# Template 28 (Easy) — Area C3: Surface-Water Hydrology
def template_rational_method_peak_flow():
    """
    Peak Discharge by the Rational Method

    Scenario:
        A small urban catchment's peak discharge follows the rational
        formula in SI form,

            Q = C * i * A / 360      (Q in m^3/s, i in mm/h, A in ha)

        The catchment is either homogeneous (single C) or composed of two
        land uses, in which case an area-weighted composite C is computed
        first — an input-structure branch.

    Difficulty: Easy
    Grounding: FHWA HEC-22, 3rd ed., Ch. 3 (rational method and runoff
        coefficients — Table 3-1, on-disk; C sampled inside each land
        use's tabulated range and stated in the question per the
        self-containment policy).
    Physical bounds: intensity i in [25, 90] mm/h with a per-sample floor
        from the Q band; single-catchment area in [max(2, 0.52/C), 20] ha,
        composite sub-areas in [3, 12] ha each; C sampled within the
        HEC-22 range of the named use; peak flow asserted in
        [0.1, 6.0] m^3/s.

    Returns:
        tuple: (question, solution)
    """
    form = random.choice(["single", "composite"])

    if form == "single":
        use = random.choice(list(_RATIONAL_USES.keys()))
        phrase = _RATIONAL_USES[use]
        c_lo, c_hi = RATIONAL_C[use]
        C = round(random.uniform(c_lo, c_hi), 2)
        A = round(random.uniform(max(2.0, 0.52 / C), 20.0), 1)
        i = random.randint(max(25, math.ceil(43.2 / (C * A))), 90)
        Q = round(C * i * A / 360.0, 2)
        setup_text = (
            f"a {A:.1f} ha catchment consisting of a {phrase} with a "
            f"runoff coefficient of C = {C:.2f} (per HEC-22)")
        steps = (
            f"**Step 1:** State the rational formula in SI form.\n"
            f"Q = C * i * A / 360, with Q in m^3/s, i in mm/h, and A in "
            f"ha.\n\n"
            f"**Step 2:** Substitute the given values.\n"
            f"Q = {C:.2f} * {i} * {A:.1f} / 360\n\n"
            f"**Step 3:** Evaluate the peak discharge.\n"
            f"Q = {Q:.2f} m^3/s\n\n"
        )
    else:
        use1, use2 = random.sample(list(_RATIONAL_USES.keys()), 2)
        p1, p2 = _RATIONAL_USES[use1], _RATIONAL_USES[use2]
        c1 = round(random.uniform(*RATIONAL_C[use1]), 2)
        c2 = round(random.uniform(*RATIONAL_C[use2]), 2)
        A1 = round(random.uniform(3.0, 12.0), 1)
        A2 = round(random.uniform(3.0, 12.0), 1)
        A = round(A1 + A2, 1)
        Cw = round((c1 * A1 + c2 * A2) / A, 3)
        i = random.randint(max(25, math.ceil(43.2 / (Cw * A))), 90)
        Q = round(Cw * i * A / 360.0, 2)
        C = Cw
        setup_text = (
            f"a catchment with two land uses: {A1:.1f} ha of {p1} "
            f"(C = {c1:.2f}) and {A2:.1f} ha of {p2} (C = {c2:.2f}), the "
            f"runoff coefficients taken from HEC-22 Table 3-1")
        steps = (
            f"**Step 1:** Compute the area-weighted composite runoff "
            f"coefficient.\n"
            f"C_w = (C1*A1 + C2*A2) / (A1 + A2) "
            f"= ({c1:.2f} * {A1:.1f} + {c2:.2f} * {A2:.1f}) / {A:.1f} "
            f"= {Cw:.3f}\n\n"
            f"**Step 2:** State the rational formula in SI form.\n"
            f"Q = C_w * i * A / 360, with Q in m^3/s, i in mm/h, and A "
            f"in ha.\n\n"
            f"**Step 3:** Substitute and evaluate.\n"
            f"Q = {Cw:.3f} * {i} * {A:.1f} / 360 = {Q:.2f} m^3/s\n\n"
        )

    assert 0.1 <= Q <= 6.0, f"peak flow out of bounds: {Q}"

    question = (
        f"During a design storm with rainfall intensity i = {i} mm/h, "
        f"determine the peak discharge in m^3/s from {setup_text}, using "
        f"the rational method in SI form (Q = C*i*A/360 with i in mm/h "
        f"and A in hectares)."
    )

    solution = (
        f"**Given:**\n"
        f"Rainfall intensity (i): {i} mm/h\n"
        f"Catchment: {setup_text}\n\n"
        f"{steps}"
        f"**Answer:** The peak discharge is {Q:.2f} m^3/s"
    )

    return question, solution


# Template 29 (Intermediate) — Area C3: Surface-Water Hydrology
def template_scs_curve_number_runoff():
    """
    Direct Runoff Depth by the SCS Curve-Number Method (SI Storm Data)

    Scenario:
        The TR-55 runoff equation is UNIT-BOUND: it is valid only with
        all depths in inches. The storm data arrive in millimeters, so
        the chain requires deliberate unit handling around the empirical
        core:

            P_in = P_mm / 25.4
            S = 1000/CN - 10        (inches)
            Ia = 0.2*S; validity requires P_in > Ia
            Q_in = (P_in - 0.2*S)^2 / (P_in + 0.8*S)
            Q_mm = Q_in * 25.4

    Difficulty: Intermediate
    Grounding: NRCS TR-55 (June 1986), Ch. 2, Eqs. 2-2/2-3/2-4 and
        Tables 2-2a/2-2c (on-disk; CN stated in the question per the
        self-containment policy; the equation's inch-bound validity is
        stated, making the unit conversions part of the required
        reasoning).
    Physical bounds: (land use, HSG, CN) from a fixed list of TR-55 rows
        with CN >= 61; storm depth sampled (in mm) inside a per-sample
        window solved from the runoff equation so the runoff lands in
        ~[0.28, 5.3] inches; runoff asserted in [0.2, 5.5] in
        (equivalently [5, 140] mm).

    Returns:
        tuple: (question, solution)
    """
    use, hsg, CN = random.choice(_SCS_COMBOS)
    assert SCS_CURVE_NUMBERS[use][{"A": 0, "B": 1, "C": 2, "D": 3}[hsg]] == CN, (
        "combo out of sync with constants table")
    phrase = _SCS_DISPLAY[use]

    S = round(1000.0 / CN - 10.0, 3)
    Ia = round(0.2 * S, 3)

    def _p_for(Qt):
        bq = 0.4 * S + Qt
        cq = 0.04 * S ** 2 - 0.8 * S * Qt
        return (bq + math.sqrt(bq ** 2 - 4 * cq)) / 2.0

    P_lo_in = max(2.0, _p_for(0.28))
    P_hi_in = min(6.5, _p_for(5.3))
    assert P_lo_in < P_hi_in, f"empty P window: {P_lo_in}, {P_hi_in}"
    P_mm = int(round(random.uniform(P_lo_in, P_hi_in) * 25.4))

    # Round-then-recompute from the presented storm depth in mm.
    P_in = round(P_mm / 25.4, 2)
    Q_in = round((P_in - Ia) ** 2 / (P_in + 0.8 * S), 3)
    Q_mm = round(Q_in * 25.4, 1)

    assert P_in > Ia, f"storm below initial abstraction: {P_in} vs {Ia}"
    assert 0.2 <= Q_in <= 5.5, f"runoff depth out of bounds: {Q_in}"

    question = (
        f"A catchment consists of {phrase} on hydrologic soil group "
        f"{hsg}, for which TR-55 gives a curve number CN = {CN}. A "
        f"design storm delivers P = {P_mm} mm of rainfall. The SCS "
        f"runoff equation is valid only with all depths in inches "
        f"(1 in = 25.4 mm): S = 1000/CN - 10 and "
        f"Q = (P - 0.2S)^2 / (P + 0.8S), applicable when P exceeds "
        f"the initial abstraction Ia = 0.2S. Determine the direct "
        f"runoff depth in millimeters."
    )

    solution = (
        f"**Given:**\n"
        f"Land use: {phrase}; HSG {hsg}; CN = {CN}\n"
        f"Storm depth (P): {P_mm} mm\n\n"
        f"**Step 1:** Convert the storm depth to inches (the SCS "
        f"equation is inch-bound).\n"
        f"P = {P_mm} / 25.4 = {P_in:.2f} in\n\n"
        f"**Step 2:** Compute the potential maximum retention.\n"
        f"S = 1000/CN - 10 = 1000/{CN} - 10 = {S:.3f} in\n\n"
        f"**Step 3:** Check the initial abstraction.\n"
        f"Ia = 0.2*S = 0.2 * {S:.3f} = {Ia:.3f} in; since P = "
        f"{P_in:.2f} in > Ia, runoff occurs and the runoff equation "
        f"applies.\n\n"
        f"**Step 4:** Compute the direct runoff depth in inches.\n"
        f"Q = (P - 0.2S)^2 / (P + 0.8S) "
        f"= ({P_in:.2f} - {Ia:.3f})^2 / ({P_in:.2f} + 0.8 * {S:.3f}) "
        f"= {Q_in:.3f} in\n\n"
        f"**Step 5:** Convert the runoff back to millimeters.\n"
        f"Q = {Q_in:.3f} * 25.4 = {Q_mm:.1f} mm\n\n"
        f"**Answer:** The direct runoff depth is {Q_mm:.1f} mm"
    )

    return question, solution


# Template 30 (Advanced) — Area C3: Surface-Water Hydrology
def template_linear_reservoir_routing_step():
    """
    Two Routing Steps Through a Linear Detention Reservoir

    Scenario:
        A detention basin's storage is proportional to its outflow,
        S = k*O. Combining this with level-pool continuity over one time
        step and solving for the unknown end-of-step outflow gives the
        working equation

            O_next = [ (I_j + I_j+1)*dt/2 + O_j*(k - dt/2) ] / (k + dt/2)

        The trace constructs this equation from continuity and the
        storage relation (showing the expansion and collection), then
        applies it over TWO consecutive intervals of a rising inflow
        hydrograph to obtain the outflow at the end of the second
        interval.

    Difficulty: Advanced
    Grounding: Chow, Maidment & Mays, Applied Hydrology, Ch. 8
        (lumped/level-pool routing typology; linear-reservoir storage
        relation) — mechanics fully self-contained in the question.
    Physical bounds: storage constant k in [1800, 5400] s with
        k >= dt/2 + 300 s; time step dt in {900, 1200, 1800} s; rising
        inflows I1 in [2, 10], I2 = I1 + [2, 10], I3 = I2 + [1.5, 8]
        m^3/s; initial outflow O1 in [0.3, 0.6]*I1; asserted
        O1 < O2 < O3 < I3 (attenuated rising outflow).

    Returns:
        tuple: (question, solution)
    """
    dt = random.choice([900, 1200, 1800])
    k = random.randint(max(1800, dt // 2 + 300) // 60, 5400 // 60) * 60
    I1 = round(random.uniform(2.0, 10.0), 1)
    I2 = round(I1 + random.uniform(2.0, 10.0), 1)
    I3 = round(I2 + random.uniform(1.5, 8.0), 1)
    O1 = round(random.uniform(0.3, 0.6) * I1, 1)

    half_dt = round(dt / 2.0, 1)
    den = round(k + half_dt, 1)
    coeff = round(k - half_dt, 1)

    # Interval 1 -> O2, then interval 2 -> O3, each from displayed values.
    num1a = round((I1 + I2) * half_dt, 1)
    num1b = round(O1 * coeff, 1)
    O2 = round((num1a + num1b) / den, 2)
    num2a = round((I2 + I3) * half_dt, 1)
    num2b = round(O2 * coeff, 1)
    O3 = round((num2a + num2b) / den, 2)

    assert coeff > 0, f"storage constant too small: {k}, {dt}"
    assert O1 < O2 < O3 < I3, (
        f"attenuation violated: {O1}, {O2}, {O3}, {I3}")

    question = (
        f"A detention basin behaves as a linear reservoir whose storage "
        f"is proportional to its outflow, S = k*O with k = {k} s. An "
        f"inflow hydrograph rises through I1 = {I1:.1f}, I2 = {I2:.1f}, "
        f"and I3 = {I3:.1f} m^3/s at successive times separated by "
        f"dt = {dt} s, and the outflow at the first time is O1 = "
        f"{O1:.1f} m^3/s. Starting from the level-pool continuity "
        f"equation (I_j + I_{{j+1}})/2 * dt - (O_j + O_{{j+1}})/2 * dt = "
        f"S_{{j+1}} - S_j, derive the working equation for the "
        f"end-of-interval outflow and apply it over both intervals to "
        f"determine the outflow O3, in m^3/s, at the end of the second "
        f"interval."
    )

    solution = (
        f"**Given:**\n"
        f"Storage relation: S = k*O, k = {k} s; dt = {dt} s\n"
        f"Inflows: I1 = {I1:.1f}, I2 = {I2:.1f}, I3 = {I3:.1f} m^3/s; "
        f"initial outflow O1 = {O1:.1f} m^3/s\n\n"
        f"**Step 1:** Substitute the storage relation into continuity "
        f"and expand.\n"
        f"(I_j + I_{{j+1}})/2 * dt - (O_j + O_{{j+1}})/2 * dt = "
        f"k*(O_{{j+1}} - O_j)\n"
        f"(I_j + I_{{j+1}})*dt/2 - O_j*dt/2 - O_{{j+1}}*dt/2 = "
        f"k*O_{{j+1}} - k*O_j\n\n"
        f"**Step 2:** Collect the unknown O_{{j+1}} on one side.\n"
        f"(I_j + I_{{j+1}})*dt/2 + O_j*(k - dt/2) = O_{{j+1}}*(k + dt/2)\n"
        f"O_{{j+1}} = [ (I_j + I_{{j+1}})*dt/2 + O_j*(k - dt/2) ] / "
        f"(k + dt/2)\n"
        f"Here k - dt/2 = {k} - {half_dt:.1f} = {coeff:.1f} s and "
        f"k + dt/2 = {k} + {half_dt:.1f} = {den:.1f} s.\n\n"
        f"**Step 3:** Evaluate the interval-1 terms (I1 -> I2).\n"
        f"(I1 + I2)*dt/2 = ({I1:.1f} + {I2:.1f}) * {half_dt:.1f} "
        f"= {num1a:.1f} m^3\n"
        f"O1*(k - dt/2) = {O1:.1f} * {coeff:.1f} = {num1b:.1f} m^3\n\n"
        f"**Step 4:** Compute the outflow at the end of interval 1.\n"
        f"O2 = ({num1a:.1f} + {num1b:.1f}) / {den:.1f} = {O2:.2f} "
        f"m^3/s\n\n"
        f"**Step 5:** Evaluate the interval-2 terms (I2 -> I3), using "
        f"the O2 just computed.\n"
        f"(I2 + I3)*dt/2 = ({I2:.1f} + {I3:.1f}) * {half_dt:.1f} "
        f"= {num2a:.1f} m^3\n"
        f"O2*(k - dt/2) = {O2:.2f} * {coeff:.1f} = {num2b:.1f} m^3\n\n"
        f"**Step 6:** Compute the outflow at the end of interval 2.\n"
        f"O3 = ({num2a:.1f} + {num2b:.1f}) / {den:.1f} = {O3:.2f} "
        f"m^3/s\n\n"
        f"**Step 7:** Confirm the routed outflow behaves physically.\n"
        f"O1 = {O1:.1f} < O2 = {O2:.2f} < O3 = {O3:.2f} m^3/s, all below "
        f"the peak inflow I3 = {I3:.1f} m^3/s — a rising, attenuated "
        f"outflow as expected for a detention basin.\n\n"
        f"**Answer:** The outflow at the end of the second interval is "
        f"{O3:.2f} m^3/s"
    )

    return question, solution
