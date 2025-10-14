import random
import numpy as np
import math
from data.templates.branches.chemical_engineering.constants import GAS_PHASE_REACTANTS, THERMO_SUBSTANCES, CRITICAL_PROPERTIES


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

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the gas volume.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    gas_name = random.choice(GAS_PHASE_REACTANTS)
    # Moles of gas
    n = round(random.uniform(0.5, 5.0), 2)
    # Temperature in Kelvin
    T_k = round(random.uniform(273.15, 500.0), 2)
    # Pressure in Pascals (for calculation)
    P_pa = round(random.uniform(100000, 500000))
    # Pressure in Kilopascals (for the question text)
    P_kpa = round(P_pa / 1000, 1)

    # Define the ideal gas constant in SI units
    R = 8.314  # Pa·m³/(mol·K)

    # 2. Perform the core calculation
    # Volume will be in cubic meters (m³)
    V_m3 = (n * R * T_k) / P_pa
    # Convert volume to Liters for the final answer
    V_L = V_m3 * 1000

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
    # Specific volume of saturated liquid in m³/kg
    V_l = round(random.uniform(0.001, 0.002), 5)
    # Specific volume of saturated vapor in m³/kg
    V_v = round(random.uniform(0.05, 2.0), 3)
    # Quality (mass fraction of vapor)
    x = round(random.uniform(0.1, 0.9), 2)

    # 2. Perform the core calculation
    V = (1 - x) * V_l + x * V_v

    # 3. Generate the question and solution strings
    question = (
        f"A closed vessel contains a saturated mixture of {substance} at a constant pressure. "
        f"The specific volume of the saturated liquid is {V_l} m³/kg and the "
        f"specific volume of the saturated vapor is {V_v} m³/kg. If the quality "
        f"of the mixture is {x*100:.0f}%, what is the overall specific volume of the mixture?"
    )

    solution = (
        f"**Step 1:** State the formula for the specific volume of a saturated mixture.\n"
        f"The formula is V = (1-x)V^l + xV^v, where V^l is the saturated liquid specific volume and V^v is the saturated vapor specific volume.\n\n"

        f"**Step 2:** List the given values.\n"
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

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the compressibility factor.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs
    R = 0.08314  # L·bar/(mol·K)

    substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
    properties = CRITICAL_PROPERTIES[substance_name]
    Tc = properties["Tc"]
    Pc = properties["Pc"]
    omega = properties["omega"]

    # Generate random T and P in the gas phase (Tr > 1) and at low-moderate pressure
    Tr = round(random.uniform(1.1, 3.0), 3)
    Pr = round(random.uniform(0.1, 2.0), 3)
    T = round(Tr * Tc, 2)
    P = round(Pr * Pc, 2)

    # 2. Perform the core calculation
    # Calculate virial coefficients
    B0 = 0.083 - 0.422 / (Tr**1.6)
    B1 = 0.139 - 0.172 / (Tr**4.2)
    # Calculate compressibility factor
    Z = 1 + (Pr / Tr) * (B0 + omega * B1)
    # Optional extension: Calculate molar volume
    V = (Z * R * T) / P

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
        f"Tr = T / Tc = {T} K / {Tc} K = {round(Tr, 4)}\n"
        f"Pr = P / Pc = {P} bar / {Pc} bar = {round(Pr, 4)}\n\n"

        f"**Step 3:** Calculate the virial equation coefficients, B0 and B1.\n"
        f"B0 = 0.083 - (0.422 / {round(Tr, 4)}**1.6) = {round(B0, 4)}\n"
        f"B1 = 0.139 - (0.172 / {round(Tr, 4)}**4.2) = {round(B1, 4)}\n\n"

        f"**Step 4:** Substitute all values to calculate the compressibility factor (Z).\n"
        f"Z = 1 + ({round(Pr, 4)} / {round(Tr, 4)}) * ({round(B0, 4)} + {omega} * {round(B1, 4)})\n"
        f"Z = 1 + {round(Pr / Tr, 4)} * ({round(B0 + omega * B1, 4)})\n"
        f"Z = {round(Z, 4)}\n\n"

        f"(Optional) **Step 5:** Calculate the molar volume (V).\n"
        f"V = Z * R * T / P = ({round(Z, 4)} * {R} * {T}) / {P} = {round(V, 4)} L/mol\n\n"

        f"**Answer:** The compressibility factor, Z, for {substance_name} at the given conditions is **{round(Z, 4)}**."
    )

    return question, solution


# Template 6 (Advanced)
def template_vdw_solve_for_volume():
    """
    Molar Volume from the Van der Waals Equation

    Scenario:
        Solving for molar volume (V) from a cubic equation of state, given
        temperature (T) and pressure (P), requires finding the roots of a
        cubic polynomial, which is done with a numerical solver.

        When the state (T, P) is below the critical point, the equation can yield
        three positive real roots, corresponding to the saturated liquid volume, the
        saturated vapor volume, and an unstable intermediate root.

        The governing equation is cast into a polynomial form for root-finding:
            V**3 + c2*V**2 + c1*V + c0 = 0
        Where:
        - c2 = -(b + R*T/P)
        - c1 = a/P
        - c0 = -(a*b)/P

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the possible molar volumes.
            - str: A step-by-step solution showing the calculation and root analysis.
    """
    # 1. Parameterize the inputs
    R = 0.08314  # L·bar/(mol·K)

    substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
    properties = CRITICAL_PROPERTIES[substance_name]
    Tc = properties["Tc"]
    Pc = properties["Pc"]

    # Choose a T and P below the critical point
    Tr = random.uniform(0.85, 0.95)
    T = round(Tr * Tc, 2)
    Pr = random.uniform(Tr * 0.8, Tr * 0.9) # Estimate a realistic saturation pressure
    P = round(Pr * Pc, 2)

    # 2. Perform the core calculation
    a = (27 * (R**2) * (Tc**2)) / (64 * Pc)
    b = (R * Tc) / (8 * Pc)

    # Coefficients of the cubic polynomial: V³ + c₂V² + c₁V + c₀ = 0
    c2 = -(b + R * T / P)
    c1 = a / P
    c0 = -(a * b) / P
    coeffs = [1, c2, c1, c0]

    roots = np.roots(coeffs)

    # Improved Root Handling: Filter for nearly-real, positive roots
    tolerance = 1e-9
    positive_real_roots = sorted([
        root.real for root in roots if abs(root.imag) < tolerance and root.real > 0
    ])

    V_f, V_g = 0, 0
    if len(positive_real_roots) >= 2: # Typically 3 for subcritical
        V_f = positive_real_roots[0]
        V_g = positive_real_roots[-1] # Safely choose the largest
    elif len(positive_real_roots) == 1:
        V_g = positive_real_roots[0]

    # 3. Generate the question and solution strings
    question = (
        f"A vessel of {substance_name} is held at a temperature of {T} K and a pressure of {P} bar. "
        f"Using the van der Waals equation of state, determine the possible molar volumes (L/mol) "
        f"for the liquid and/or vapor phases. The critical properties for {substance_name} are:\n"
        f"Tc = {Tc} K\nPc = {Pc} bar"
    )

    solution = (
        f"**Step 1:** Calculate the van der Waals parameters 'a' and 'b'.\n"
        f"a = (27 * R² * Tc²) / (64 * Pc) = {round(a, 4)} L²·bar/mol²\n"
        f"b = (R * Tc) / (8 * Pc) = {round(b, 5)} L/mol\n\n"

        f"**Step 2:** Formulate the cubic equation: V³ + c₂V² + c₁V + c₀ = 0.\n"
        f"The calculated coefficients, with their respective units, are:\n"
        f"- c₂ = {round(c2, 4)} (L/mol)\n"
        f"- c₁ = {round(c1, 4)} (L²/mol²)\n"
        f"- c₀ = {round(c0, 6)} (L³/mol³)\n\n"

        f"**Step 3:** Solve the polynomial for its roots using a numerical solver.\n"
        f"The physically meaningful (positive, real) roots found are: "
        f"{', '.join([f'{r:.4f}' for r in positive_real_roots])} L/mol\n\n"

        f"**Step 4:** Interpret the physical significance of the roots.\n"
    )

    if len(positive_real_roots) >= 2:
        solution += (
            f"For a subcritical state, we expect multiple positive real roots.\n"
            f"- The smallest root corresponds to the molar volume of the **liquid phase (Vf)**.\n"
            f"- The largest root corresponds to the molar volume of the **vapor phase (Vg)**.\n"
            f"- Any intermediate root lies on the thermodynamically unstable branch of the isotherm (where pressure incorrectly increases with volume) and is disregarded.\n\n"
            f"**Answer:**\n"
            f"  - Saturated Liquid Volume (Vf) ≈ **{round(V_f, 4)} L/mol**\n"
            f"  - Saturated Vapor Volume (Vg) ≈ **{round(V_g, 4)} L/mol**"
        )
    elif len(positive_real_roots) == 1:
        solution += (
            f"A single positive real root indicates the substance exists in a single phase (gas or supercritical fluid).\n\n"
            f"**Answer:**\n"
            f"  - Molar Volume (V) = **{round(V_g, 4)} L/mol**"
        )
    else:
        solution += "No physically meaningful (positive, real) roots were found for these conditions, which may indicate an issue with the applicability of the model at this state."

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

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the work of compression.
            - str: A step-by-step solution showing the calculation and comparison
                   to the ideal gas case.
    """
    # 1. Parameterize the inputs
    R = 0.08314  # L·bar/(mol·K)

    substance_name = random.choice(list(CRITICAL_PROPERTIES.keys()))
    properties = CRITICAL_PROPERTIES[substance_name]
    Tc = properties["Tc"]
    Pc = properties["Pc"]
    omega = properties["omega"]

    # Generate conditions in the gas phase (Tr > 1) at moderate pressures
    Tr = round(random.uniform(1.2, 3.0), 3)
    T = round(Tr * Tc, 2)
    
    P1_r = round(random.uniform(0.5, 2.0), 2)
    P2_r = round(random.uniform(2.5, 5.0), 2)
    P1 = round(P1_r * Pc, 2)
    P2 = round(P2_r * Pc, 2)

    # 2. Perform the core calculation
    # Calculate B
    B0 = 0.083 - 0.422 / (Tr**1.6)
    B1 = 0.139 - 0.172 / (Tr**4.2)
    B = (R * Tc / Pc) * (B0 + omega * B1) # Units: L/mol

    # Calculate initial and final states
    Z1 = 1 + (B * P1) / (R * T)
    Z2 = 1 + (B * P2) / (R * T)
    V1 = (Z1 * R * T) / P1
    V2 = (Z2 * R * T) / P2

    # Calculate work using the integrated virial equation
    term1 = R * T * math.log(V2 / V1)
    term2 = -B * R * T * ((1/V2) - (1/V1))
    W_virial_Lbar = -(term1 + term2)
    W_virial_J = W_virial_Lbar * 100  # Convert L·bar to Joules

    # For comparison, calculate ideal gas work
    W_ideal_Lbar = -R * T * math.log(P1 / P2)
    W_ideal_J = W_ideal_Lbar * 100

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
        f"The integrated form is: W = -[RT·ln(V2/V1) - BRT(1/V2 - 1/V1)]\n\n"

        f"**Step 2:** Calculate the second virial coefficient (B) at T = {T} K.\n"
        f"Reduced Temperature, Tr = T/Tc = {T}/{Tc} = {Tr}\n"
        f"B0 = 0.083 - 0.422 / ({Tr})**1.6 = {round(B0, 4)}\n"
        f"B1 = 0.139 - 0.172 / ({Tr})**4.2 = {round(B1, 4)}\n"
        f"B = (R·Tc/Pc) * (B0 + ω·B1) = {round(B, 5)} L/mol\n\n"

        f"**Step 3:** Determine the initial (V1) and final (V2) molar volumes.\n"
        f"Z1 = 1 + B·P1/(R·T) = 1 + ({round(B, 5)}*{P1})/({R}*{T}) = {round(Z1, 4)}\n"
        f"V1 = Z1·R·T/P1 = {round(V1, 5)} L/mol\n"
        f"Z2 = 1 + B·P2/(R·T) = 1 + ({round(B, 5)}*{P2})/({R}*{T}) = {round(Z2, 4)}\n"
        f"V2 = Z2·R·T/P2 = {round(V2, 5)} L/mol\n\n"

        f"**Step 4:** Substitute V1 and V2 into the integrated work equation.\n"
        f"W = -[{round(R*T, 2)}·ln({round(V2, 5)}/{round(V1, 5)}) - {round(B*R*T, 3)}(1/{round(V2, 5)} - 1/{round(V1, 5)})]\n"
        f"W = -[{round(term1, 2)} + {round(term2, 2)}] = {round(W_virial_Lbar, 2)} L·bar/mol\n\n"

        f"**Step 5:** Convert the work to the required units (J/mol).\n"
        f"Since 1 L·bar = 100 J:\n"
        f"W = {round(W_virial_Lbar, 2)} L·bar/mol * 100 J/(L·bar) = {round(W_virial_J, 0)} J/mol\n\n"

        f"**For Comparison:** The work required for an ideal gas is W_ideal = -RT·ln(P1/P2) = {round(W_ideal_J, 0)} J/mol. "
        f"The deviation shows the effect of intermolecular forces accounted for by the virial equation.\n\n"

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
