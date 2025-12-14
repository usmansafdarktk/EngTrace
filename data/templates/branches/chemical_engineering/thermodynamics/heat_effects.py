import random
from scipy.optimize import fsolve
from data.templates.branches.chemical_engineering.constants import SUBSTANCES_FOR_HEATING, SUBSTANCES_FOR_VAPORIZATION, HEATS_OF_FORMATION, REACTIONS, CP_PARAMS, COMBUSTION_REACTIONS


# Template 1 (Easy)
def template_sensible_heat_constant_cp():
    """
    Sensible Heat Calculation (Constant Heat Capacity)

    Scenario:
        This template calculates the sensible heat required to change the
        temperature of a substance without changing its phase. It uses the
        substance's actual specific heat capacity (Cp) for a physically
        accurate problem.

        The governing equation on a mass basis is:
            Q = m * Cp * (T2 - T1)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the sensible heat.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs using the substance list
    substance_data = random.choice(SUBSTANCES_FOR_HEATING)
    substance_name = substance_data["name"]
    substance_state = substance_data["state"]
    Cp = substance_data["Cp"]  # J/(g·K)
    
    # Retrieve safe temperature limits from the constant data
    # Defaults provided in case keys are missing in legacy data
    min_temp = substance_data.get("min_temp", 20)
    max_temp = substance_data.get("max_temp", 150)

    # Generate a random mass in grams
    m = round(random.uniform(50.0, 1000.0), 1)

    # Generate temperatures within the safe bounds for this specific substance
    # Ensure there's at least a 10 degree window for T2
    safe_max_T1 = max(min_temp, max_temp - 15) 
    T1_C = round(random.uniform(min_temp, safe_max_T1), 1)
    
    # Ensure T2 is higher than T1 but within max limit
    T2_C = round(random.uniform(T1_C + 10, max_temp), 1)

    # 2. Perform the core calculation
    delta_T = T2_C - T1_C
    # Q will be in Joules (g * J/g·K * K)
    Q_J = m * Cp * delta_T
    # Convert to kilojoules for the final answer
    Q_kJ = Q_J / 1000

    # 3. Generate the question and solution strings
    question = (
        f"Calculate the heat in kJ required to raise the temperature of {m} g "
        f"of {substance_state} {substance_name} from {T1_C}°C to {T2_C}°C. "
        f"The specific heat capacity of {substance_name} is {Cp} J/g·K. "
        f"Assume constant pressure and no phase change occurs."
    )

    solution = (
        f"**Step 1:** State the formula.\n"
        f"Q = m * Cp * ΔT\n\n"
        f"**Note:** We assume the specific heat capacity (Cp) is constant over this temperature range.\n"
        f"\n\n"

        f"**Step 2:** List the given values.\n"
        f"- Mass (m) = {m} g\n"
        f"- Specific Heat Capacity (Cp) = {Cp} J/g·K\n"
        f"- Initial Temperature (T1) = {T1_C}°C\n"
        f"- Final Temperature (T2) = {T2_C}°C\n\n"

        f"**Step 3:** Calculate the temperature difference (ΔT).\n"
        f"A change in temperature has the same magnitude in Celsius and Kelvin.\n"
        f"ΔT = T2 - T1 = {T2_C}°C - {T1_C}°C = {round(delta_T, 1)} K\n\n"

        f"**Step 4:** Substitute the values into the formula to find the heat in Joules (J).\n"
        f"Q = {m} g * {Cp} J/g·K * {round(delta_T, 1)} K\n"
        f"Q = {round(Q_J, 1)} J\n\n"

        f"**Step 5:** Convert the heat to kilojoules (kJ) as requested.\n"
        f"Q = {round(Q_J, 1)} J * (1 kJ / 1000 J) = {round(Q_kJ, 2)} kJ\n\n"

        f"**Answer:** The heat required is **{round(Q_kJ, 2)} kJ**."
    )

    return question, solution


# Template 2 (Easy)
def template_latent_heat_vaporization():
    """
    Latent Heat of Vaporization

    Scenario:
        This template calculates the heat required to cause a phase change
        (vaporization) at a constant temperature and pressure. This energy,
        known as latent heat, is used to overcome intermolecular forces
        rather than to increase the substance's kinetic energy (temperature).

        The governing equation is:
            Q = n * ΔH_vap
        Where:
        - Q: Total heat absorbed
        - n: Number of moles
        - ΔH_vap: Molar heat of vaporization

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the latent heat.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs using the substance list
    substance_data = random.choice(SUBSTANCES_FOR_VAPORIZATION)
    substance_name = substance_data["name"]
    delta_H_vap = substance_data["delta_H_vap"]  # in kJ/mol

    # Generate a random number of moles
    n = round(random.uniform(0.5, 5.0), 2)  # moles

    # 2. Perform the core calculation
    # The result will be in kJ since ΔH_vap is in kJ/mol
    Q_kJ = n * delta_H_vap

    # 3. Generate the question and solution strings
    question = (
        f"How much heat in kJ is required to completely vaporize {n} moles of "
        f"liquid {substance_name} at its normal boiling point? The molar heat of "
        f"vaporization for {substance_name} is {delta_H_vap} kJ/mol."
    )

    solution = (
        f"**Step 1:** State the formula.\n"
        f"Q = n * ΔH_vap\n\n"

        f"**Step 2:** List the given values.\n"
        f"- Moles (n) = {n} mol\n"
        f"- Molar Heat of Vaporization (ΔH_vap) = {delta_H_vap} kJ/mol\n\n"

        f"**Step 3:** Substitute the values into the formula.\n"
        f"Q = {n} mol * {delta_H_vap} kJ/mol\n\n"

        f"**Step 4:** Calculate the total heat required.\n"
        f"Q = {round(Q_kJ, 2)} kJ\n\n"

        f"**Answer:** The total heat required to vaporize the {substance_name} is **{round(Q_kJ, 2)} kJ**."
    )

    return question, solution


# Template 3 (Easy)
def template_heat_of_reaction_formation():
    """
    Standard Heat of Reaction from Heats of Formation

    Scenario:
        This template applies Hess's Law to find the standard heat of reaction
        at 298.15 K (ΔH_rxn°). It is calculated by subtracting the sum of the
        standard heats of formation (ΔH_f°) of the reactants from the sum of
        the heats of formation of the products, with each being weighted by
        its stoichiometric coefficient (v).

        The governing equation is:
            ΔH_rxn° = Σ(v * ΔH_f°)products - Σ(v * ΔH_f°)reactants

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the standard heat of reaction.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs by choosing a random reaction
    reaction_data = random.choice(REACTIONS)
    equation = reaction_data["equation"]
    reactants = reaction_data["reactants"]
    products = reaction_data["products"]

    # 2. Perform the core calculation with full precision
    products_enthalpy = sum(nu * HEATS_OF_FORMATION[species] for species, nu in products.items())
    reactants_enthalpy = sum(nu * HEATS_OF_FORMATION[species] for species, nu in reactants.items())
        
    # Final calculation uses the precise, unrounded intermediate values
    delta_H_rxn = products_enthalpy - reactants_enthalpy

    # 3. Generate the question and solution strings
    required_species = list(reactants.keys()) + list(products.keys())
    hf_list_for_question = "\n".join(
        f"- {s}: {HEATS_OF_FORMATION[s]} kJ/mol" for s in sorted(list(set(required_species)))
    )

    # Question updated for clarity
    question = (
        f"Calculate the standard enthalpy of reaction, ΔH_rxn°, at 298.15 K (in kJ) for the reaction as written:\n"
        f"{equation}\n\n"
        f"Use the standard heats of formation (ΔH_f°) provided below:\n{hf_list_for_question}"
    )

    # Helper functions to format the solution steps
    def format_sum(species_dict):
        return " + ".join([f"({nu} * ΔH_f°[{s}])" for s, nu in species_dict.items()])

    def format_calc(species_dict):
        return " + ".join([f"({nu} * {HEATS_OF_FORMATION[s]})" for s, nu in species_dict.items()])

    # Solution updated to use correct rounding procedures and units
    solution = (
        f"**Step 1:** State the formula.\n"
        f"ΔH_rxn° = Σ(v * ΔH_f°)products - Σ(v * ΔH_f°)reactants\n\n"

        f"**Step 2:** Calculate the total enthalpy of the products.\n"
        f"Σ_products = {format_sum(products)}\n"
        f"Σ_products = {format_calc(products)}\n"
        f"Σ_products = {products_enthalpy:.3f} kJ\n\n"

        f"**Step 3:** Calculate the total enthalpy of the reactants.\n"
        f"Σ_reactants = {format_sum(reactants)}\n"
        f"Σ_reactants = {format_calc(reactants)}\n"
        f"Σ_reactants = {reactants_enthalpy:.3f} kJ\n\n"

        f"**Step 4:** Calculate the standard enthalpy of reaction using the full-precision values.\n"
        f"ΔH_rxn° = ({products_enthalpy:.3f}) - ({reactants_enthalpy:.3f})\n"
        f"ΔH_rxn° = {round(delta_H_rxn, 2)} kJ\n\n"

        f"**Answer:** The standard enthalpy of reaction is **{round(delta_H_rxn, 2)} kJ** for the reaction as written."
    )
    
    return question, solution


# Template 4 (Intermediate)
def template_sensible_heat_temp_dependent_cp():
    """
    Sensible Heat (Temperature-Dependent Heat Capacity)

    Scenario:
        This template provides a more accurate calculation of sensible heat where
        the heat capacity (Cp) is a polynomial function of temperature. The total
        heat required is found by integrating this function over the temperature
        range from T1 to T2.

        The governing equation is:
            Q = n * integral(Cp(T) dT) from T1 to T2

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the sensible heat via integration.
            - str: A step-by-step solution showing the analytical integration.
    """
    # 1. Parameterize the inputs
    R = 8.314  # J/(mol·K)
    substance_name = random.choice(list(CP_PARAMS.keys()))
    params = CP_PARAMS[substance_name]
    A, B, C, D = params["A"], params["B"], params["C"], params["D"]

    n = round(random.uniform(1.0, 10.0), 2)  # moles

    # Determine valid temperature ranges based on phase to ensure physical plausibility
    if "(l)" in substance_name:
        # Liquids: Keep range lower to avoid boiling (approximate general cap)
        T1 = round(random.uniform(280.0, 300.0), 2)
        T2 = round(random.uniform(T1 + 20, 350.0), 2) 
    else:
        # Gases/Solids: Can handle higher temperatures
        T1 = round(random.uniform(298.15, 500.0), 2)
        T2 = round(random.uniform(T1 + 100, 1200.0), 2)

    # 2. Perform the core calculation via analytical integration
    # integral(A + BT + CT² + DT⁻²)dT = AT + (B/2)T² + (C/3)T³ - D/T
    def integral_mean_cp_over_r(T):
        term_a = A * T
        term_b = (B / 2) * (T**2)
        term_c = (C / 3) * (T**3)
        term_d = -D / T if T != 0 else 0
        return term_a + term_b + term_c + term_d

    # Calculate the definite integral per mole
    val_T2 = integral_mean_cp_over_r(T2)
    val_T1 = integral_mean_cp_over_r(T1)
    integral_val = R * (val_T2 - val_T1)
    
    Q_J = n * integral_val  # Total heat in Joules
    Q_kJ = Q_J / 1000       # Convert to kilojoules

    # 3. Generate the question and solution strings
    # Dynamically create the Cp/R equation string for the question
    cp_eq_parts = [str(A)]
    if B != 0: cp_eq_parts.append(f"{B:.3e}*T")
    if C != 0: cp_eq_parts.append(f"{C:.3e}*T^2")
    if D != 0: cp_eq_parts.append(f"{D:.3e}*T^-2")
    
    cp_eq_str = " + ".join(cp_eq_parts).replace("+ -", "- ")

    question = (
        f"Calculate the heat in kJ required to raise the temperature of {n} moles of "
        f"{substance_name} from {T1} K to {T2} K. The molar heat capacity for "
        f"{substance_name} is given by the equation:\n"
        f"  Cp/R = {cp_eq_str}\n"
        f"Where the constants are: A={A}, B={B:.3g}, C={C:.3g}, D={D:.3g}"
    )

    # Helper to format the polynomial substitution string for the solution
    def format_poly_sub(T_val):
        terms = [f"{A}*({T_val})"]
        if B != 0: terms.append(f"({B:.3e}/2)*({T_val})^2")
        if C != 0: terms.append(f"({C:.3e}/3)*({T_val})^3")
        if D != 0: terms.append(f"-({D:.3e}/{T_val})")
        return " + ".join(terms).replace("+ -", "- ")

    solution = (
        f"**Step 1:** State the integral formula for total heat (Q).\n"
        f"Q = n * R * ∫[from {T1} to {T2}] (Cp/R) dT\n\n"

        f"**Step 2:** Perform the analytical integration.\n"
        f"The indefinite integral of Cp(T)/R is:\n"
        f"I(T) = A*T + (B/2)*T^2 + (C/3)*T^3 - D/T\n\n"

        f"**Step 3:** Evaluate the definite integral.\n"
        f"I({T2}) = {format_poly_sub(T2)} = {val_T2:.4f}\n"
        f"I({T1}) = {format_poly_sub(T1)} = {val_T1:.4f}\n"
        f"Integral Value = I({T2}) - I({T1}) = {val_T2:.4f} - {val_T1:.4f} = {val_T2 - val_T1:.4f} K\n\n"

        f"**Step 4:** Calculate the total heat Q.\n"
        f"Q = n * R * (Integral Value)\n"
        f"Q = {n} mol * 8.314 J/(mol K) * {val_T2 - val_T1:.4f} K\n"
        f"Q = {Q_J:.2f} J\n\n"

        f"**Step 5:** Convert to kilojoules.\n"
        f"Q = {Q_J:.2f} J / 1000 = {Q_kJ:.2f} kJ\n\n"

        f"**Answer:** The heat required is **{Q_kJ:.2f} kJ**."
    )

    return question, solution


# Template 5 (Advanced)
def template_adiabatic_flame_temperature():
    """
    Adiabatic Flame Temperature

    Scenario:
        This template calculates the theoretical maximum temperature the products
        of a combustion reaction can reach if it occurs adiabatically (no heat loss).
        This requires an energy balance: the heat released by the reaction at a
        reference temperature (298.15 K) must be completely absorbed by the
        product gases as sensible heat, raising their temperature from the
        reference temperature to the final adiabatic flame temperature (T_ad).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the adiabatic flame temperature.
            - str: A step-by-step solution showing the numerical method.
    """
    # 1. Parameterize inputs
    R = 8.314  # J/(mol·K)
    T_initial = 298.15 # K
    
    # Validation Loop for Physics Sanity Check
    while True:
        try:
            reaction_data = random.choice(COMBUSTION_REACTIONS)
            fuel = reaction_data["fuel"]
            equation = reaction_data["equation"]
            reactants = reaction_data["reactants"]
            products = reaction_data["products"]

            # 2. Perform the core calculation
            # Calculate the standard heat of reaction at 298.15 K
            products_enthalpy_298 = sum(nu * HEATS_OF_FORMATION[s] for s, nu in products.items())
            reactants_enthalpy_298 = sum(nu * HEATS_OF_FORMATION[s] for s, nu in reactants.items())
            delta_H_298_kJ = products_enthalpy_298 - reactants_enthalpy_298
            delta_H_298_J = delta_H_298_kJ * 1000 # Convert kJ to J

            # Define integral of Cp/R
            def integral_mean_cp_over_r(T, A, B, C, D):
                return A*T + (B/2)*T**2 + (C/3)*T**3 - D/T

            # Define the energy balance equation (Sensible Heat + Heat of Rxn = 0)
            def energy_balance(T_final):
                sensible_heat_products = 0
                for species, nu in products.items():
                    params = CP_PARAMS[species]
                    A, B, C, D = params["A"], params["B"], params["C"], params["D"]
                    
                    # Integral from T_initial to T_final
                    val_T = integral_mean_cp_over_r(T_final, A, B, C, D)
                    val_T0 = integral_mean_cp_over_r(T_initial, A, B, C, D)
                    
                    sensible_heat_products += nu * R * (val_T - val_T0)
                    
                return sensible_heat_products + delta_H_298_J

            # Solve for the root
            initial_guess = 2000 
            adiabatic_temp_kelvin = fsolve(energy_balance, initial_guess)[0]

            #  PHYSICS SANITY CHECK 
            # Reject results that are physically implausible (< 1500 K or > 3500 K)
            # This filters out bad data or unit mismatches that cause "broken" problems.
            if 1500 < adiabatic_temp_kelvin < 3500:
                break # Valid result found, exit loop

        except (KeyError, ValueError, IndexError):
            # If data is missing (KeyError) or solver fails, retry with a different reaction
            continue

    # 3. Generate the question and solution strings
    question = (
        f"{fuel} gas enters a furnace at {T_initial} K and is burned completely with "
        f"the theoretical amount of dry air (also at {T_initial} K). Assuming the "
        f"process is adiabatic and there is no shaft work, estimate the "
        f"adiabatic flame temperature in Kelvin."
    )

    solution = (
        f"**Step 1:** Write the balanced chemical equation including inert nitrogen from air.\n"
        f"  {equation}\n\n"

        f"**Step 2:** Calculate the standard heat of reaction at {T_initial} K (ΔH_rxn°).\n"
        f"Using standard heats of formation:\n"
        f"ΔH_rxn° = {delta_H_298_kJ:.2f} kJ\n\n"

        f"**Step 3:** Set up the energy balance equation.\n"
        f"For an adiabatic process, the heat released by the reaction must be absorbed as sensible heat by the products.\n"
        f"Heat Absorbed by Products + Heat of Reaction = 0\n"
        f"  Σ [ n_i * ∫(Cp_i dT) from {T_initial} to T_ad ] + ΔH_rxn° = 0\n\n"

        f"**Step 4:** Solve the equation for the adiabatic flame temperature (T_ad).\n"
        f"This equation is non-linear because Cp is a function of temperature. We use a numerical solver to find T_ad.\n"
        f"Initial guess: {initial_guess} K\n\n"

        f"**Step 5:** State the result.\n"
        f"The calculated adiabatic flame temperature is:\n"
        f"T_ad = {adiabatic_temp_kelvin:.2f} K\n\n"

        f"**Answer:** The estimated adiabatic flame temperature is **{adiabatic_temp_kelvin:.2f} K**."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each heat effects template with different random seeds
    and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/thermodynamics/heat_effects.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_sensible_heat_constant_cp, "sensible_heat_constant_cp", "Easy"),
        (template_latent_heat_vaporization, "latent_heat_vaporization", "Easy"),
        (template_heat_of_reaction_formation, "heat_of_reaction_formation", "Easy"),
        (template_sensible_heat_temp_dependent_cp, "sensible_heat_temp_dependent_cp", "Intermediate"),
        (template_adiabatic_flame_temperature, "adiabatic_flame_temperature", "Advanced"),
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
                "area": "heat_effects",
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
