import random
import numpy as np
import math
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches.chemical_engineering.constants import GAS_PHASE_REACTANTS, THERMO_SUBSTANCES, CRITICAL_PROPERTIES, REAL_FLUID_DATA
from data.templates.branches._emission import signed_term, joined_terms, paren_neg


def _hu(x, places):
    """Round half-up to `places` dp, resolving the tie in DECIMAL.

    `round()` resolves a half-way tie on the binary value and disagrees with a
    reader doing decimal arithmetic (spec P2 as amended, DECISIONS D-012).
    """
    q = Decimal(1).scaleb(-places)
    d = x if isinstance(x, Decimal) else Decimal(repr(x))
    v = d.quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


def _as_printed(x, spec):
    """The value a reader recovers from `x` when it is printed with `spec`.

    P2 asks that the stored value and the printed value be the SAME value.
    Rounding alone does not achieve that: a float one ulp away from its own
    printed form puts the template and the reader on opposite sides of a
    display tie (D-016 part 2).
    """
    return float(format(x, spec))


def _is_display_tie(x, places, rel_band=1e-12):
    """Is `x` at, or within a hair of, a half-way tie at `places` dp?

    A tie is the one case where no rounding convention is defensible - a
    decimal reader applying half-up and a binary reader applying `round()`
    disagree, and the printed line closes for only one of them. Such instances
    are resampled rather than resolved (D-016).

    A narrow BAND is quarantined rather than a point, because two independent
    evaluations of the same exact quantity differ by a few ulps and an exact
    rational tie lands on opposite sides of them.
    """
    scaled = abs(x) * 10.0 ** places
    band = max(scaled * rel_band, 1e-9)
    return abs((scaled - math.floor(scaled)) - 0.5) <= band


def _decimals(x):
    """Decimal places of the shortest decimal that round-trips `x`.

    0.1 -> 1, 0.012 -> 3, 4.003 -> 3, 32.0 -> 0. Used to print a quantity
    built from table values at the precision it actually has (D-037), never
    coarser.
    """
    exponent = Decimal(repr(float(x))).normalize().as_tuple().exponent
    return max(0, -exponent)


# Template 1 (Easy)
def template_ideal_gas_volume():
    """
    Ideal Gas Law Volume Calculation

    Scenario:
        The Ideal Gas Law, PV = nRT, is a fundamental equation of state that
        describes the behavior of many gases. In this scenario, the pressure,
        temperature, and number of moles of a gas are provided. The objective
        is to apply the Ideal Gas Law to calculate the volume the gas occupies
        using the formula:

            V = nRT / P

    Screen pass 1 (2026-09-23):
        One judge found Step 2 equating the stated 1-dp kPa pressure to an
        integer Pa value from an independent rounding (`179.2 kPa = 179214
        Pa`; 99% of draws). The Pa value is now the stated kPa times 1000
        exactly and the volume is computed from it (D-016 part 2); the
        question text is unchanged and the gold volume moves in its last
        digit on 32% of seeds. V in m^3 is bound through its 5-dp display
        before the litre conversion, and a draw whose 5-dp display sits on
        a half-way tie is redrawn rather than resolved (D-016).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the gas volume.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    gas_name = random.choice(GAS_PHASE_REACTANTS)
    # Define the ideal gas constant in SI units
    R = 8.314  # Pa·m³/(mol·K)

    for _attempt in range(200):
        # Moles of gas
        n = round(random.uniform(0.5, 5.0), 2)
        # Temperature in Kelvin
        T_k = round(random.uniform(273.15, 500.0), 2)
        # Pressure: sampled in Pa, STATED in the question in kPa at 1 dp. The
        # Pa value the solution computes with is the stated value times 1000
        # exactly, not the sampled integer it was rounded from (D-016 part 2).
        P_kpa = round(round(random.uniform(100000, 500000)) / 1000, 1)
        P_pa = int(round(P_kpa * 1000))

        # 2. Perform the core calculation
        # Volume will be in cubic meters (m³); it is printed at 5 dp and then
        # converted, so it is bound through that display, and a draw whose
        # 5-dp display sits on a half-way tie is redrawn (D-016).
        V_m3 = (n * R * T_k) / P_pa
        if _is_display_tie(V_m3, 5):
            continue
        V_m3 = _as_printed(V_m3, '.5f')
        # Convert volume to Liters for the final answer
        V_L = V_m3 * 1000
        break
    else:
        raise RuntimeError(
            "ideal_gas_volume: no display-stable sample in 200 draws")

    # 3. Generate the question and solution strings
    question = (
        f"Calculate the volume in liters occupied by {n} moles of {gas_name} gas "
        f"at a temperature of {T_k} K and a pressure of {P_kpa} kPa. "
        f"Assume the gas behaves ideally."
    )

    solution = (
        f"**Step 1:** State the Ideal Gas Law formula solved for volume (V).\n"
        f"The formula is V = \\frac{{nRT}}{{P}}\n\n"

        f"**Step 2:** List the given values and the ideal gas constant, ensuring consistent SI units.\n"
        f"- Moles (n) = {n} mol\n"
        f"- Temperature (T) = {T_k} K\n"
        f"- Pressure (P) = {P_kpa} kPa = {P_pa} Pa\n"
        f"- Ideal Gas Constant (R) = {R} Pa·m³/(mol·K)\n\n"

        f"**Step 3:** Substitute the values into the equation.\n"
        f"V = \\frac{{({n} \\text{{ mol}}) \\times ({R} \\text{{ Pa·m³/mol·K}}) \\times ({T_k} \\text{{ K}})}}{{{P_pa} \\text{{ Pa}}}}\n\n"

        f"**Step 4:** Calculate the volume in cubic meters.\n"
        f"V = {round(V_m3, 5)} \\text{{ m³}}\n\n"

        f"**Step 5:** Convert the volume to liters as requested.\n"
        f"Since 1 \\text{{ m³}} = 1000 \\text{{ L}}:\n"
        f"V = {round(V_m3, 5)} \\text{{ m³}} \\times 1000 \\frac{{\\text{{L}}}}{{\\text{{m³}}}} = {round(V_L, 2)} \\text{{ L}}\n\n"

        f"**Answer:** The volume occupied by the gas is **{round(V_L, 2)} liters**."
    )

    return question, solution


# Template 2 (Easy)
def template_two_phase_specific_volume():
    """
    Specific Volume of a Two-Phase Mixture

    Scenario:
        When a substance exists as a mixture of liquid and vapor in equilibrium
        (a saturated mixture), its overall specific volume (V) depends on the
        proportions of the liquid and vapor phases. This proportion is defined
        by the quality (x), which is the mass fraction of vapor. This template
        provides the specific volumes of the saturated liquid (V^l) and
        saturated vapor (V^v) along with the quality, and asks to calculate
        the overall specific volume.

        The governing equation is:
            V = (1-x)V^l + xV^v

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the mixture's specific volume.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    substance = random.choice(THERMO_SUBSTANCES)
    # C3.3: the substance's own saturated volumes, at the temperature
    # REAL_FLUID_DATA states them for. They were invented with random.uniform,
    # so a named fluid could be given any volumes at all. Removing those two
    # draws moves the quality drawn next on every seed (P6).
    sat = REAL_FLUID_DATA[substance]
    T_C = sat["temp_C"]
    # Specific volume of saturated liquid in m³/kg
    V_l = sat["v_f"]
    # Specific volume of saturated vapor in m³/kg
    V_v = sat["v_g"]
    # Quality (mass fraction of vapor)
    x = round(random.uniform(0.1, 0.9), 2)

    # 2. Perform the core calculation
    V = (1 - x) * V_l + x * V_v

    # 3. Generate the question and solution strings
    question = (
        f"A closed vessel contains a saturated mixture of {substance} at {T_C} °C. "
        f"The specific volume of the saturated liquid is {V_l} m³/kg and the "
        f"specific volume of the saturated vapor is {V_v} m³/kg. If the quality "
        f"of the mixture is {x*100:.0f}%, what is the overall specific volume of the mixture?"
    )

    solution = (
        f"**Step 1:** State the formula for the specific volume of a saturated mixture.\n"
        f"The formula is V = (1-x)V^l + xV^v, where V^l is the saturated liquid specific volume and V^v is the saturated vapor specific volume.\n\n"

        f"**Step 2:** List the given values.\n"
        f"- Temperature = {T_C} °C (saturation)\n"
        f"- Quality (x) = {x}\n"
        f"- Saturated liquid specific volume (V^l) = {V_l} m³/kg\n"
        f"- Saturated vapor specific volume (V^v) = {V_v} m³/kg\n\n"

        f"**Step 3:** Substitute the values into the equation.\n"
        f"V = (1 - {x})({V_l} \\text{{ m³/kg}}) + ({x})({V_v} \\text{{ m³/kg}})\n"
        f"V = ({(1-x):.2f})({V_l}) + ({x})({V_v})\n"
        f"V = {round((1-x)*V_l, 5)} + {round(x*V_v, 5)}\n\n"

        f"**Step 4:** Calculate the final specific volume.\n"
        f"V = {round(V, 5)} \\text{{ m³/kg}}\n\n"

        f"**Answer:** The overall specific volume of the mixture is **{round(V, 5)} m³/kg**."
    )

    return question, solution


# Template 3 (Easy)
def template_rackett_equation_volume():
    """
    Rackett Equation for Saturated Liquid Volume

    Scenario:
        The Rackett equation is a generalized correlation used to estimate the
        molar volume of a saturated liquid when experimental data is not
        available. It relies on the substance's critical properties.

        This template provides the critical temperature (Tc), critical volume (Vc),
        and critical compressibility factor (Zc) for a substance. The goal is to
        calculate the saturated liquid molar volume (V_sat) at a given
        temperature (T).

        The governing equation is:
            V_sat = Vc * Zc**((1 - Tr)**0.2857)
        Where:
        - V_sat: Molar volume of the saturated liquid (cm³/mol)
        - Vc: Critical volume (cm³/mol)
        - Zc: Critical compressibility factor
        - Tr: Reduced temperature (T / Tc)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the saturated liquid volume.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
    properties = CRITICAL_PROPERTIES[substance_name]
    Tc = properties["Tc"]
    Vc = properties["Vc"]
    Zc = properties["Zc"]

    # Generate a random temperature below the critical temperature
    T = round(random.uniform(0.5 * Tc, 0.95 * Tc), 2)

    # 2. Perform the core calculation
    # Calculate the reduced temperature
    Tr = T / Tc
    # Apply the Rackett equation
    exponent = (1 - Tr)**0.2857
    V_sat = Vc * (Zc**exponent)

    # 3. Generate the question and solution strings
    question = (
        f"Estimate the molar volume of saturated liquid {substance_name} at {T} K "
        f"using the Rackett equation. The critical properties for {substance_name} are:\n"
        f"Tc = {Tc} K\n"
        f"Vc = {Vc} cm³/mol\n"
        f"Zc = {Zc}"
    )

    solution = (
        f"**Step 1:** State the Rackett equation.\n"
        f"V_sat = Vc * Zc**((1 - Tr)**0.2857)\n\n"

        f"**Step 2:** List the given properties and calculate the reduced temperature (Tr).\n"
        f"- Critical Temperature (Tc) = {Tc} K\n"
        f"- Critical Volume (Vc) = {Vc} cm³/mol\n"
        f"- Critical Compressibility (Zc) = {Zc}\n"
        f"- Temperature (T) = {T} K\n\n"
        f"Tr = T / Tc = {T} / {Tc} = {round(Tr, 4)}\n\n"

        f"**Step 3:** Substitute the values into the Rackett equation.\n"
        f"V_sat = {Vc} * {Zc}**((1 - {round(Tr, 4)})**0.2857)\n"
        f"V_sat = {Vc} * {Zc}**({round(1 - Tr, 4)}**0.2857)\n"
        f"V_sat = {Vc} * {Zc}**({round(exponent, 4)})\n"
        f"V_sat = {Vc} * {round(Zc**exponent, 4)}\n\n"

        f"**Step 4:** Calculate the final molar volume.\n"
        f"V_sat = {round(V_sat, 2)} cm³/mol\n\n"

        f"**Answer:** The estimated molar volume of saturated liquid {substance_name} at {T} K is **{round(V_sat, 2)} cm³/mol**."
    )

    return question, solution


# Template 4 (Intermediate)
def template_vdw_solve_for_pressure():
    """
    Van der Waals Equation for Pressure Calculation

    Scenario:
        The van der Waals equation is an equation of state that improves upon the
        ideal gas law by including terms for intermolecular attraction ('a') and
        molecular volume ('b'). This makes it applicable to real fluids in both
        liquid and gas phases.

        This template provides the critical properties (Tc and Pc) for a substance,
        along with a given temperature (T) and molar volume (V). The objective is
        to first calculate the van der Waals parameters 'a' and 'b', and then
        use them to find the pressure (P).

        The governing equations are:
            P = (R * T) / (V - b) - a / (V**2)
            a = (27 * R**2 * Tc**2) / (64 * Pc)
            b = (R * Tc) / (8 * Pc)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the pressure.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs
    # Gas constant in L·bar/(mol·K)
    R = 0.08314

    substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
    properties = CRITICAL_PROPERTIES[substance_name]
    Tc = properties["Tc"]
    Pc = properties["Pc"]

    # Generate a random temperature, avoiding the critical region (0.95Tc to 1.05Tc)
    while True:
        T = round(random.uniform(0.8 * Tc, 2.5 * Tc), 2)
        # Avoid temperatures too close to critical point where VdW equation is less accurate
        if not (0.95 * Tc <= T <= 1.05 * Tc):
            break

    # First, calculate 'b' to ensure V > b
    b = (R * Tc) / (8 * Pc)
    # Generate a random molar volume greater than b
    V = round(random.uniform(1.5 * b, 150 * b), 4)

    # 2. Perform the core calculation
    # Calculate parameter 'a'
    a = (27 * (R**2) * (Tc**2)) / (64 * Pc)

    # Calculate pressure using the Van der Waals equation
    P = (R * T) / (V - b) - a / (V**2)

    # 3. Generate the question and solution strings
    question = (
        f"Using the van der Waals equation of state, calculate the pressure in bar "
        f"exerted by {substance_name} at a temperature of {T} K and a molar volume "
        f"of {V} L/mol. The critical constants for {substance_name} are:\n"
        f"Tc = {Tc} K\n"
        f"Pc = {Pc} bar"
    )

    solution = (
        f"**Step 1:** State the necessary formulas.\n"
        f"Pressure: P = (R * T) / (V - b) - a / (V**2)\n"
        f"Parameter 'a': a = (27 * R**2 * Tc**2) / (64 * Pc)\n"
        f"Parameter 'b': b = (R * Tc) / (8 * Pc)\n\n"

        f"**Step 2:** List the given values and the gas constant.\n"
        f"- Temperature (T) = {T} K\n"
        f"- Molar Volume (V) = {V} L/mol\n"
        f"- Critical Temperature (Tc) = {Tc} K\n"
        f"- Critical Pressure (Pc) = {Pc} bar\n"
        f"- Gas Constant (R) = {R} L·bar/(mol·K)\n\n"

        f"**Step 3:** Calculate the co-volume parameter 'b'.\n"
        f"b = ({R} * {Tc}) / (8 * {Pc}) = {round(b, 5)} L/mol\n\n"

        f"**Step 4:** Calculate the attraction parameter 'a'.\n"
        f"a = (27 * ({R})**2 * ({Tc})**2) / (64 * {Pc}) = {round(a, 4)} L²·bar/mol²\n\n"

        f"**Step 5:** Substitute all values into the van der Waals equation to find P.\n"
        f"P = ({R} * {T}) / ({V} - {round(b, 5)}) - {round(a, 4)} / ({V})**2\n"
        f"P = {round(R * T, 2)} / {round(V - b, 5)} - {round(a, 4)} / {round(V**2, 5)}\n"
        f"P = {round((R * T) / (V - b), 2)} - {round(a / (V**2), 2)}\n\n"

        f"**Step 6:** Calculate the final pressure.\n"
        f"P = {round(P, 2)} bar\n\n"

        f"**Answer:** The pressure exerted by the {substance_name} is **{round(P, 2)} bar**."
    )

    return question, solution


# Templat 5 (Intermediate)
def template_pitzer_correlation_z():
    """
    Compressibility Factor from Pitzer's Correlation

    Scenario:
        The Pitzer correlation is a widely used application of the principle of
        corresponding states to estimate the properties of real fluids. It
        improves upon two-parameter correlations by introducing the acentric
        factor (omega), which accounts for the non-sphericity of molecules.

        This template uses a simplified form of the Pitzer correlation, valid
        for gases at low to moderate pressures, to calculate the compressibility
        factor (Z).

        The governing equations are:
            Z = 1 + (Pr / Tr) * (B0 + omega * B1)
            B0 = 0.083 - (0.422 / Tr**1.6)
            B1 = 0.139 - (0.172 / Tr**4.2)
        Where:
        - Z: Compressibility factor
        - Tr, Pr: Reduced temperature and pressure
        - omega: Acentric factor
        - B0, B1: Second virial coefficients

    Trace integrity (Layer 0, 2026-09-23):
        The sampler still CHOOSES reduced conditions to place the instance in
        the gas phase at low-moderate pressure, but the question states T and
        P, so Tr and Pr are recomputed FORWARD from the stated values at the
        4 dp the solution prints them with; the sampled pre-image used to be
        printed (Tr = 1.986 where 10.33 K / 5.2 K gives 1.9865). Every
        intermediate that is printed and then consumed - B0, B1, Pr/Tr and Z
        at 4 dp - is bound through its display first, so the molar volume is
        computed from the stated Z (D-016 part 2). B0 + omega*B1 is exact at
        4 + dec(omega) dp for every draw and is printed at that length, so it
        never rounds (D-037, D-044); Pr/Tr is a quotient whose 4-dp display
        can sit exactly on a tie only for a terminating Tr, and such a draw is
        redrawn (D-016, D-045).

    Screen pass 1 (2026-09-23):
        One judge flagged states outside the correlation's accepted range
        (Pr up to 2 at Tr near 1.1). The validity criterion applied is Smith,
        Van Ness & Abbott's (7th ed., sec. 3.6, Fig. 3.14): the two-term
        virial equation with the generalized B-correlation is adequate where
        Vr = V/Vc >= 2, the region below the line of that figure. A draw
        outside it - 2.7% of draws over 5,000 seeds, all at Tr <= 1.52 and
        high Pr (>= 1.4) - is redrawn; the sampling ranges are unchanged,
        so only those seeds change (3.0% of 500). The high-T end of the
        range is kept: thermal stability is
        not a condition of the correlation. A second judge's claim that the
        printed B0 + omega*B1 (0.0003192) differs from the value used is
        rejected: -0.0399 + 0.304*0.1323 = 0.0003192 exactly, which is the
        7-dp display the code binds.

    Layer 2 review (2026-09-25):
        One expert of three (the other two accepted) found hydrocarbons
        sampled at 889-1064 K (propylene, n-hexane, n-pentane on three of
        the five hand-check seeds), where the real substance would crack,
        and asked for a per-substance temperature cap. The equations,
        constants and arithmetic were confirmed correct. Not adopted, for
        the reason the screen-pass paragraph above already gives: the item
        is a corresponding-states correlation exercise at a stated (T, P),
        and the correlation's validity condition (Vr >= 2) is enforced;
        the chemical lifetime of the species at that state is not a
        condition of the correlation. A cap was also simulated: with
        organics capped at 750 K (an onset-of-cracking figure that has no
        on-disk source) 42% of currently valid draws would be rejected and
        nine organics (n-hexane to p-xylene, methanol, ethanol, acetone)
        would fall below a fifth of their present share of the pool, since
        the cap collides with the Vr >= 2 filter at high Pr; at 900 K, 30%
        rejected. No number, text or sample changes.

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the compressibility factor.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs
    R = 0.08314  # L·bar/(mol·K)

    for _attempt in range(200):
        substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
        properties = CRITICAL_PROPERTIES[substance_name]
        Tc = properties["Tc"]
        Pc = properties["Pc"]
        Vc = properties["Vc"]  # cm^3/mol
        omega = properties["omega"]

        # Generate random T and P in the gas phase (Tr > 1) and at low-moderate
        # pressure. The question states T and P, so Tr and Pr are recomputed
        # FORWARD from the stated values at the 4 dp the solution prints them
        # with; the sampled pre-image no longer reaches the solution (it used
        # to: `Tr = 1.986` was printed where 10.33 K / 5.2 K gives 1.9865).
        Tr_sampled = round(random.uniform(1.1, 3.0), 3)
        Pr_sampled = round(random.uniform(0.1, 2.0), 3)
        T = round(Tr_sampled * Tc, 2)
        P = round(Pr_sampled * Pc, 2)
        Tr = _as_printed(T / Tc, '.4f')
        Pr = _as_printed(P / Pc, '.4f')

        # 2. Perform the core calculation. Every intermediate that is printed
        # and then consumed is bound through its own display first, so the
        # next line is computed from the operands the reader sees (D-016).
        B0 = _as_printed(0.083 - 0.422 / (Tr**1.6), '.4f')
        B1 = _as_printed(0.139 - 0.172 / (Tr**4.2), '.4f')
        # Pr/Tr is a quotient; its 4-dp display sits exactly on a half-way tie
        # only when Tr terminates (T = 2.0000 Tc, say). Such a draw is repeated
        # rather than resolved (D-016, D-045).
        if _is_display_tie(Pr / Tr, 4):
            continue
        ratio = _as_printed(Pr / Tr, '.4f')
        # B0 + omega*B1 is EXACT at 4 + dec(omega) dp for every draw (B0 and
        # B1 are 4-dp values, omega a table value), so it is printed at that
        # length and nothing rounds: at 4 dp it sat on a tie for 10% of ethane
        # draws (omega = 0.100) and 25% of ammonia (0.250) (D-037, D-044).
        bsum_dp = 4 + _decimals(omega)
        bsum = _as_printed(B0 + omega * B1, f'.{bsum_dp}f')
        # Calculate compressibility factor
        Z = _as_printed(1 + ratio * bsum, '.4f')
        # Optional extension: Calculate molar volume, from the STATED Z
        V = (Z * R * T) / P
        # Validity: the two-term virial equation with the Pitzer B-correlation
        # is adequate only where Vr = V/Vc >= 2, the region below the line of
        # SVA Fig. 3.14 (Smith, Van Ness & Abbott, 7th ed., sec. 3.6). A draw
        # outside it (2.7% of draws: Tr <= 1.52 with high Pr) is redrawn; the
        # sampling ranges are unchanged.
        if V * 1000.0 / Vc < 2.0:
            continue
        break
    else:
        raise RuntimeError(
            "pitzer_correlation_z: no display-stable sample in 200 draws")

    # 3. Generate the question and solution strings
    question = (
        f"For {substance_name} at a temperature of {T} K and a pressure of {P} bar, "
        f"determine the compressibility factor, Z, using the Pitzer correlation for the "
        f"second virial coefficient. The properties for {substance_name} are:\n"
        f"Critical Temperature (Tc) = {Tc} K\n"
        f"Critical Pressure (Pc) = {Pc} bar\n"
        f"Acentric Factor (ω) = {omega}"
    )

    solution = (
        f"**Step 1:** State the governing formulas.\n"
        f"Z = 1 + (Pr / Tr) * (B0 + omega * B1)\n"
        f"B0 = 0.083 - (0.422 / Tr**1.6)\n"
        f"B1 = 0.139 - (0.172 / Tr**4.2)\n\n"

        f"**Step 2:** Calculate the reduced temperature (Tr) and reduced pressure (Pr).\n"
        f"Tr = T / Tc = {T} K / {Tc} K = {Tr}\n"
        f"Pr = P / Pc = {P} bar / {Pc} bar = {Pr}\n\n"

        f"**Step 3:** Calculate the virial equation coefficients, B0 and B1.\n"
        f"B0 = 0.083 - (0.422 / {Tr}**1.6) = {B0}\n"
        f"B1 = 0.139 - (0.172 / {Tr}**4.2) = {B1}\n\n"

        f"**Step 4:** Substitute all values to calculate the compressibility factor (Z).\n"
        f"Z = 1 + ({Pr} / {Tr}) * ({B0} + {paren_neg(omega)} * {B1})\n"
        f"Z = 1 + {ratio} * ({bsum:.{bsum_dp}f})\n"
        f"Z = {Z}\n\n"

        f"(Optional) **Step 5:** Calculate the molar volume (V).\n"
        f"V = Z * R * T / P = ({Z} * {R} * {T}) / {P} = {round(V, 4)} L/mol\n\n"

        f"**Answer:** The compressibility factor, Z, for {substance_name} at the given conditions is **{Z}**."
    )

    return question, solution


# Template 6 (Advanced)
def template_vdw_solve_for_volume():
    """
    Molar Volume from the Van der Waals Equation

    Scenario:
        Solving for molar volume (V) from a cubic equation of state, given
        temperature (T) and pressure (P), requires finding the roots of a
        cubic polynomial.

        When the state (T, P) is within the two-phase region (subcritical), 
        the equation yields three real roots:
        1. Smallest root: Liquid-phase molar volume.
        2. Largest root: Vapor-phase molar volume.
        3. Intermediate root: Unstable state (physically meaningless).

        The governing equation is cast into a polynomial form for root-finding:
            V^3 - (b + RT/P)V^2 + (a/P)V - (ab/P) = 0

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the possible molar volumes.
            - str: A step-by-step solution showing the calculation and root analysis.
    """
    # 1. Parameterize the inputs with a Validation Loop
    # We want to ensure we generate a problem with 3 roots (the "Two-Phase" scenario)
    # as this is the most pedagogically valuable case for cubic EOS problems.
    
    R = 0.08314  # L·bar/(mol·K)
    
    # Loop until we find a set of inputs that produces 3 distinct real roots
    while True:
        try:
            substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
            properties = CRITICAL_PROPERTIES[substance_name]
            Tc = properties["Tc"]
            Pc = properties["Pc"]

            # Choose T well below critical point to ensure the "loop" is wide enough
            Tr = random.uniform(0.65, 0.85)
            T = round(Tr * Tc, 2)
            
            # Estimate VdW saturation pressure to hit the 3-root region
            # Approx VdW vapor pressure: log10(Pr) ~ -3(1/Tr - 1)
            approx_Pr_sat = 10**(-3.0 * (1.0/Tr - 1.0))
            
            # Pick P close to this saturation pressure
            Pr = approx_Pr_sat * random.uniform(0.95, 1.05)
            P = round(Pr * Pc, 3) # Use 3 decimals for precision

            # 2. Perform the core calculation
            a = (27 * (R**2) * (Tc**2)) / (64 * Pc)
            b = (R * Tc) / (8 * Pc)

            # Coefficients of the cubic polynomial: V³ + c₂V² + c₁V + c₀ = 0
            c2 = -(b + R * T / P)
            c1 = a / P
            c0 = -(a * b) / P
            coeffs = [1, c2, c1, c0]

            roots = np.roots(coeffs)

            # Filter for positive, real roots
            tolerance = 1e-9
            positive_real_roots = sorted([
                root.real for root in roots if abs(root.imag) < tolerance and root.real > 0
            ])
            
            # If we found 3 roots, we have a valid problem. Break the loop.
            if len(positive_real_roots) == 3:
                V_f = positive_real_roots[0]
                V_g = positive_real_roots[-1]
                break
                
        except Exception:
            continue

    # 3. Generate the question and solution strings
    question = (
        f"A vessel of {substance_name} is held at a temperature of {T} K and a pressure of {P} bar. "
        f"Using the van der Waals equation of state, determine the possible molar volumes (L/mol) "
        f"predicted for the liquid and vapor phases.\n\n"
        f"The critical properties for {substance_name} are:\n"
        f"Tc = {Tc} K\nPc = {Pc} bar"
    )

    solution = (
        f"**Step 1:** Calculate the van der Waals parameters 'a' and 'b'.\n"
        f"a = (27 * R^2 * Tc^2) / (64 * Pc) = {round(a, 4)} L^2·bar/mol^2\n"
        f"b = (R * Tc) / (8 * Pc) = {round(b, 5)} L/mol\n\n"

        f"**Step 2:** Formulate the cubic equation: V^3 + c2*V^2 + c1*V + c0 = 0.\n"
        f"The coefficients are derived from the VdW equation:\n"
        f"c2 = -(b + RT/P) = {round(c2, 4)} L/mol\n"
        f"c1 = a/P = {round(c1, 4)} L^2/mol^2\n"
        f"c0 = -(ab)/P = {round(c0, 6)} L^3/mol^3\n\n"

        f"**Step 3:** Solve the polynomial for its roots.\n"
        f"Using a numerical solver, we find three real, positive roots:\n"
        f"Roots: {', '.join([f'{r:.4f}' for r in positive_real_roots])} L/mol\n\n"
        f"\n\n"

        f"**Step 4:** Interpret the physical significance of the roots.\n"
        f"For a state in the subcritical region, the EOS predicts three roots:\n"
        f"- The **smallest root** corresponds to the molar volume of the **Liquid Phase**.\n"
        f"- The **largest root** corresponds to the molar volume of the **Vapor Phase**.\n"
        f"- The intermediate root is physically unstable and disregarded.\n\n"

        f"**Answer:**\n"
        f"1. Liquid-phase Volume (V_liq) ≈ **{round(V_f, 4)} L/mol**\n"
        f"2. Vapor-phase Volume (V_vap) ≈ **{round(V_g, 4)} L/mol**"
    )

    return question, solution


# Template 7 (Advanced)
def template_work_isothermal_virial():
    """
    Work of Isothermal Compression for a Virial Gas

    Scenario:
        Calculating the work (W) for a mechanically reversible, isothermal
        process requires integrating the pressure (P) with respect to volume (V).
        For a non-ideal gas described by the virial equation, the P-V
        relationship is more complex than the ideal gas law, leading to a
        different result for the work of compression.

        This template calculates work by first determining the second virial
        coefficient (B) and the initial/final state properties (V1, V2), and
        then applying the analytical integral of the virial equation.

        The governing equations are:
            P = R*T * (1/V + B/V**2)
            W = - integral(P dV) from V1 to V2
            W = -[R*T*ln(V2/V1) - B*R*T*(1/V2 - 1/V1)]

    Trace integrity (Layer 0, 2026-09-23):
        The question states T, P1 and P2, so Tr is recomputed FORWARD from
        the stated T at the 3 dp the solution prints it with; the sampled
        pre-image used to be printed (helium: 14.73/5.2 = 2.8327 was printed
        as the sampled 2.832). Every intermediate that is printed and then
        consumed is bound through its own display first - B0 and B1 at 4 dp,
        B at 5 dp (Z1 = 1 + B*P1/(R*T) did not close from the printed 5-dp B),
        Z1 and Z2 at 4 dp, V1 and V2 at 5 dp, R*T at 2 dp, B*R*T at 3 dp and
        the two work terms at 2 dp - so W is the sum of the printed terms
        (D-016 part 2). R*T and B*R*T are finite decimals whose displays can
        sit exactly on a half-way tie (R*T = 20.785 at 250.00 K, nitrogen, is
        the one such point on the input grid); such a draw is redrawn (D-016).

    Screen pass 1 (2026-09-23):
        One judge found the solution mixing virial forms: Z1 and Z2 came
        from the pressure-explicit Z = 1 + B*P/(R*T), while the work
        integral of Step 1 is that of the volume-explicit P = RT(1/V +
        B/V^2). Confirmed. Under the pressure-explicit form V - RT/P = B is
        constant, so the isothermal work between two pressures is
        RT*ln(P2/P1) exactly (the ideal-gas value), and the printed
        deviation from it was an artefact of the mix (mean -7 J/mol, down
        to -104, over the valid states; down to -14,000 outside them). The
        volume-explicit form is now used throughout: Z is the physical root
        of Z^2 - Z - B*P/(R*T) = 0, V = Z*R*T/P, and Step 4 is unchanged.
        The deviation is then second order in B and positive (mean +6
        J/mol on the valid states). Every gold answer changes. The same
        judge's point about states where the two-term equation is not valid
        is taken with the SVA criterion Vr = V/Vc >= 2 at both states (7th
        ed., sec. 3.6, Fig. 3.14): about half of all draws fail it at the
        final state P2 (sampled up to 5 Pc) - over 5,000 seeds 31% by
        Vr < 2 and 19% with no real root - and are redrawn; the sampling
        ranges are unchanged. Item-pool effect: the kept pool sits at higher
        Tr (mean 2.07 -> 2.47, 10th percentile 1.38 -> 2.00) and slightly
        lower P2/Pc (mean 3.73 -> 3.52); all 27 substances remain, the light
        gases (ethane, ethylene, CO2, N2, O2) depleted by about a third.
        Another judge's claim that
        Z1 = 1.022 and V1 = 0.81278 did not follow from the printed operands
        is rejected: 0.01746*129.3/(0.08314*1236.83) = 0.021954, so
        Z1 = 1.0220, and 1.022*0.08314*1236.83/129.3 = 0.81278.

    Layer 2 review (2026-09-25):
        One expert of three (the other two accepted) found the organic
        substances compressed at 900-1188 K on all five hand-check seeds
        (n-pentane 1045 K, acetone 1188 K, n-hexane 1088 K), where the
        real species would crack during a sustained reversible process,
        and asked for a per-substance temperature cap; the expert also
        noted, correctly, that the Vr >= 2 filter skews the kept pool to
        higher Tr (mean 2.47). Equations, constants and arithmetic were
        confirmed correct. Not adopted: the item is an equation-of-state
        exercise whose stated validity condition (Vr >= 2 at both states)
        is enforced, and no per-substance thermal-stability limit exists
        on disk to source a cap from. A cap was simulated: with organics
        capped at 750 K, 57% of currently valid draws would be rejected
        and ten of the 27 substances (n-pentane through p-xylene, the
        aromatics, methanol, ethanol, acetone) would leave the pool
        entirely, because at Tr <= 1.5 no P2 in the sampled 2.5-5 Pc
        window satisfies Vr >= 2; at 900 K four substances would still be
        lost. Removing a third of the substance table on one reviewer's
        note is a sampling-range decision for the owner, not a fix; it is
        recorded here and in the Layer 2 report. No number, text or sample
        changes.

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the work of compression.
            - str: A step-by-step solution showing the calculation and comparison
                   to the ideal gas case.
    """
    # 1. Parameterize the inputs
    R = 0.08314  # L·bar/(mol·K)

    for _attempt in range(200):
        substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
        properties = CRITICAL_PROPERTIES[substance_name]
        Tc = properties["Tc"]
        Pc = properties["Pc"]
        omega = properties["omega"]
        Vc = properties["Vc"]  # cm^3/mol

        # Generate conditions in the gas phase (Tr > 1) at moderate pressures.
        # The question states T, so Tr is recomputed FORWARD from the stated T
        # at the 3 dp the solution prints it with; the sampled pre-image no
        # longer reaches the solution (helium: 14.73/5.2 = 2.8327 was printed
        # as the sampled 2.832).
        Tr_sampled = round(random.uniform(1.2, 3.0), 3)
        T = round(Tr_sampled * Tc, 2)

        P1_r = round(random.uniform(0.5, 2.0), 2)
        P2_r = round(random.uniform(2.5, 5.0), 2)
        P1 = round(P1_r * Pc, 2)
        P2 = round(P2_r * Pc, 2)
        Tr = _as_printed(T / Tc, '.3f')

        # 2. Perform the core calculation. Every intermediate that is printed
        # and then consumed is bound through its own display first (D-016):
        # the next line is computed from the operands the reader sees.
        B0 = _as_printed(0.083 - 0.422 / (Tr**1.6), '.4f')
        B1 = _as_printed(0.139 - 0.172 / (Tr**4.2), '.4f')
        B = _as_printed((R * Tc / Pc) * (B0 + omega * B1), '.5f')  # Units: L/mol

        # Initial and final states from the SAME equation the work integral
        # uses, the volume-explicit Z = 1 + B/V (screen pass 1): with
        # V = Z*R*T/P it is Z^2 - Z - B*P/(R*T) = 0, whose physical root is
        # Z = (1 + sqrt(1 + 4*B*P/(R*T)))/2. Before, Z came from the
        # pressure-explicit Z = 1 + B*P/(R*T), under which V - R*T/P = B is
        # constant and the isothermal work is R*T*ln(P2/P1) exactly, so the
        # printed deviation from the ideal gas was an artefact of the mix.
        # A state with no real root is outside the validity region below.
        disc1 = 1 + 4 * B * P1 / (R * T)
        disc2 = 1 + 4 * B * P2 / (R * T)
        if disc1 <= 0 or disc2 <= 0:
            continue
        # Z is printed at 4 dp and V at 5 dp, then consumed; each is bound
        # through its display, and a draw whose display sits on a half-way
        # tie is redrawn (D-016).
        Z1_raw = (1 + math.sqrt(disc1)) / 2
        Z2_raw = (1 + math.sqrt(disc2)) / 2
        if _is_display_tie(Z1_raw, 4) or _is_display_tie(Z2_raw, 4):
            continue
        Z1 = _as_printed(Z1_raw, '.4f')
        Z2 = _as_printed(Z2_raw, '.4f')
        if _is_display_tie(Z1 * R * T / P1, 5) or _is_display_tie(Z2 * R * T / P2, 5):
            continue
        V1 = _as_printed((Z1 * R * T) / P1, '.5f')
        V2 = _as_printed((Z2 * R * T) / P2, '.5f')
        # Validity: the two-term virial equation is adequate only where
        # Vr = V/Vc >= 2, the region below the line of SVA Fig. 3.14 (Smith,
        # Van Ness & Abbott, 7th ed., sec. 3.6), at BOTH states. About half of
        # all draws fail it at P2 (sampled up to 5 Pc); they are redrawn and
        # the sampling ranges are unchanged.
        if V1 * 1000.0 / Vc < 2.0 or V2 * 1000.0 / Vc < 2.0:
            continue

        # R*T (7 dp) and B*R*T (12 dp) are finite decimals, so their 2-dp and
        # 3-dp displays CAN sit exactly on a half-way tie: R*T = 20.785 at
        # T = 250.00 K (nitrogen) is the one such point on the input grid.
        # Such a draw is repeated rather than resolved (D-016).
        if _is_display_tie(R * T, 2) or _is_display_tie(B * R * T, 3):
            continue
        RT = _as_printed(R * T, '.2f')
        BRT = _as_printed(B * R * T, '.3f')

        # Calculate work using the integrated virial equation
        term1 = _as_printed(RT * math.log(V2 / V1), '.2f')
        term2 = _as_printed(-BRT * ((1/V2) - (1/V1)), '.2f')
        W_virial_Lbar = _as_printed(-(term1 + term2), '.2f')
        W_virial_J = W_virial_Lbar * 100  # Convert L·bar to Joules

        # For comparison, calculate ideal gas work
        W_ideal_Lbar = -R * T * math.log(P1 / P2)
        W_ideal_J = W_ideal_Lbar * 100
        break
    else:
        raise RuntimeError(
            "work_isothermal_virial: no display-stable sample in 200 draws")

    # 3. Generate the question and solution strings
    question = (
        f"Calculate the work in J/mol required to isothermally and reversibly "
        f"compress 1 mole of {substance_name} from {P1} bar to {P2} bar at a "
        f"constant temperature of {T} K. Base your calculation on the virial "
        f"equation of state truncated to two terms. The properties for {substance_name} are:\n"
        f"Tc = {Tc} K, Pc = {Pc} bar, ω = {omega}"
    )

    solution = (
        f"**Step 1:** Define the pressure from the virial equation and find the analytical integral for work.\n"
        f"P = RT(1/V + B/V²)\n"
        f"The integrated form is: W = -(RT·ln(V2/V1) - BRT·(1/V2 - 1/V1))\n\n"

        f"**Step 2:** Calculate the second virial coefficient (B) at T = {T} K.\n"
        f"Reduced Temperature, Tr = T/Tc = {T}/{Tc} = {Tr}\n"
        f"B0 = 0.083 - 0.422 / ({Tr})**1.6 = {B0}\n"
        f"B1 = 0.139 - 0.172 / ({Tr})**4.2 = {B1}\n"
        f"B = (R·Tc/Pc) * (B0 + ω·B1) = {B} L/mol\n\n"

        f"**Step 3:** Determine the initial (V1) and final (V2) molar volumes from the same equation, Z = 1 + B/V. "
        f"With V = Z·R·T/P it becomes Z² - Z - B·P/(R·T) = 0, whose physical root is Z = (1 + sqrt(1 + 4·B·P/(R·T)))/2.\n"
        f"Z1 = (1 + sqrt(1 + 4·B·P1/(R·T)))/2 = (1 + sqrt(1 + 4*{paren_neg(B)}*{P1}/({R}*{T})))/2 = {Z1}\n"
        f"V1 = Z1·R·T/P1 = {Z1}*{R}*{T}/{P1} = {V1} L/mol\n"
        f"Z2 = (1 + sqrt(1 + 4·B·P2/(R·T)))/2 = (1 + sqrt(1 + 4*{paren_neg(B)}*{P2}/({R}*{T})))/2 = {Z2}\n"
        f"V2 = Z2·R·T/P2 = {Z2}*{R}*{T}/{P2} = {V2} L/mol\n\n"

        f"**Step 4:** Substitute V1 and V2 into the integrated work equation.\n"
        f"W = -({RT}·ln({V2}/{V1}) {signed_term(-BRT)}*(1/{V2} - 1/{V1}))\n"
        f"W = -({term1} {signed_term(term2)}) = {W_virial_Lbar} L·bar/mol\n\n"

        f"**Step 5:** Convert the work to the required units (J/mol).\n"
        f"Since 1 L·bar = 100 J:\n"
        f"W = {W_virial_Lbar} L·bar/mol * 100 J/(L·bar) = {round(W_virial_J, 0)} J/mol\n\n"

        f"**For Comparison:** The work required for an ideal gas is W_ideal = -RT·ln(P1/P2) = {round(W_ideal_J, 0)} J/mol. "
        f"The two-term virial equation changes the work only at second order in B: to first order V - RT/P = B is a constant offset that does no work between fixed pressures, so the deviation is small.\n\n"

        f"**Answer:** The required work of compression is approximately **{round(W_virial_J, 0)} J/mol**."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each volumetric properties of pure fluids template with 
    different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/thermodynamics/volumetric_properties_pure_fluids.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_ideal_gas_volume, "ideal_gas_volume", "Easy"),
        (template_two_phase_specific_volume, "two_phase_specific_volume", "Easy"),
        (template_rackett_equation_volume, "rackett_equation_volume", "Easy"),
        (template_vdw_solve_for_pressure, "vdw_solve_for_pressure", "Intermediate"),
        (template_pitzer_correlation_z, "pitzer_correlation_z", "Intermediate"),
        (template_vdw_solve_for_volume, "vdw_solve_for_volume", "Advanced"),
        (template_work_isothermal_virial, "work_isothermal_virial", "Advanced"),
    ]

    # List to store all generated problems
    all_problems = []

    # Generate problems for each template
    for template_func, id_name, level in templates:
        for _ in range(50):
            # Generate a unique seed for each problem
            seed = random.randint(1_000_000_000, 4_000_000_000)
            random.seed(seed)

            # Generate the problem and solution
            question, solution = template_func()

            # Create a JSON entry
            problem_entry = {
                "seed": seed,
                "branch": "chemical_engineering",
                "domain": "thermodynamics",
                "area": "volumetric_properties_pure_fluids",
                "id": id_name,
                "level": level,
                "question": question,
                "solution": solution
            }

            # Add to the list of problems
            all_problems.append(problem_entry)

    # Shuffle the problems to mix templates and levels
    random.shuffle(all_problems)

    # Write all problems to a .jsonl file
    with open(output_file, "w") as file:
        for problem in all_problems:
            file.write(json.dumps(problem))
            file.write("\n")

    print(f"Successfully generated {len(all_problems)} problems and saved to {output_file}")


if __name__ == "__main__":
    main()
