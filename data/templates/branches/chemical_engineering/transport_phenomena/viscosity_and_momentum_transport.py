import random
import math
from data.templates.branches.chemical_engineering.constants import COMMON_LIQUIDS, COMMON_GASES, GAS_MOLECULAR_PARAMS, POWER_LAW_FLUIDS


# Template 1 (Easy)
def template_newtons_law_shear_stress():
    """
    Shear Stress and Force Calculation for Flow Between Parallel Plates

    Scenario:
        This template models Couette flow, a classic fluid dynamics problem where a layer
        of fluid is sheared between two parallel plates. One plate is stationary, and
        the other moves at a constant velocity. Assuming a linear velocity profile,
        we use Newton's Law of Viscosity to find the shear stress in the fluid and
        the force required to move the plate.

        The relevant equations are:
            - Velocity Gradient: dvx/dy = V / Y
            - Shear Stress: tau_yx = mu * (dvx/dy)
            - Force: F = tau_yx * A

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute shear stress and force.
            - str: A step-by-step solution showing the calculations.
    """
    # 1. Parameterize the inputs with random values
    fluid_name, (density, viscosity) = random.choice(list(COMMON_LIQUIDS.items()))
    V = round(random.uniform(0.1, 2.0), 2)
    Y_cm = round(random.uniform(0.1, 2.5), 2)
    Y_m = Y_cm / 100
    A = round(random.uniform(0.5, 5.0), 2)

    # 2. Perform the core calculation
    velocity_gradient = V / Y_m
    shear_stress = viscosity * velocity_gradient
    force = shear_stress * A

    # 3. Generate the question and solution strings
    question = (
        f"Two large parallel plates with an area of {A} m^2 each are separated by a "
        f"thin film of {fluid_name} that is {Y_cm} cm thick. The top plate is moved "
        f"at a constant velocity of {V} m/s, while the bottom plate is held stationary.\n\n"
        f"Assuming the fluid exhibits Newtonian behavior and a linear velocity profile, calculate:\n"
        f"a) The shear stress (tau_yx) exerted on the fluid.\n"
        f"b) The total force (F) required to move the top plate at the given velocity."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Dynamic Viscosity of {fluid_name} (mu): {viscosity:.2e} Pa·s\n"
        f"- Plate Area (A): {A} m^2\n"
        f"- Plate Velocity (V): {V} m/s\n"
        f"- Distance between plates (Y): {Y_cm} cm = {Y_m} m\n\n"
        
        f"**Step 1:** Calculate the velocity gradient (dvx/dy).\n"
        f"For a linear velocity profile between a stationary and a moving plate, the gradient is constant:\n"
        f"dvx/dy = V / Y\n"
        f"dvx/dy = {V} m/s / {Y_m} m = {velocity_gradient:.2f} 1/s\n\n"
        
        f"**Step 2:** Calculate the shear stress (tau_yx).\n"
        f"Using Newton's Law of Viscosity:\n"
        f"tau_yx = mu * (dvx/dy)\n"
        f"tau_yx = ({viscosity:.2e} Pa·s) * ({velocity_gradient:.2f} 1/s)\n"
        f"tau_yx = {shear_stress:.3f} Pa\n\n"
        
        f"**Step 3:** Calculate the force (F).\n"
        f"Force is the shear stress acting over the entire area of the plate:\n"
        f"F = tau_yx * A\n"
        f"F = ({shear_stress:.3f} Pa) * ({A} m^2)\n"
        f"F = {force:.3f} N\n\n"
        
        f"**Answer:**\n"
        f"a) The shear stress in the fluid is **{shear_stress:.3f} Pa**.\n"
        f"b) The force required to move the plate is **{force:.3f} N**."
    )

    return question, solution


# Template 2 (Easy)
def template_kinematic_viscosity():
    """
    Kinematic Viscosity Calculation

    Scenario:
        This template tests the definition of kinematic viscosity (ν), which is the
        ratio of the dynamic viscosity (μ) to the density (ρ) of a fluid. It is a
        measure of a fluid's internal resistance to flow under gravitational forces.

        The relevant equations are:
            - Kinematic Viscosity: ν = μ / ρ
            - Unit Conversion: 1 m²/s = 10,000 Stokes (St)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the kinematic viscosity.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    
    # Combine liquids and gases into one dictionary for random selection
    all_fluids = {**COMMON_LIQUIDS, **COMMON_GASES}
    fluid_name, (density, viscosity) = random.choice(list(all_fluids.items()))
    
    # Randomly choose the target units for the answer
    target_units = random.choice(["m^2/s", "Stokes (St)"])

    # 2. Perform the core calculation
    
    # Calculate kinematic viscosity in SI units (m²/s)
    nu_si = viscosity / density
    
    # Perform unit conversion if necessary
    final_nu = nu_si
    if "Stokes" in target_units:
        final_nu = nu_si * 10000

    # 3. Generate the question and solution strings
    
    question = (
        f"The dynamic viscosity (μ) of {fluid_name} at standard conditions is approximately "
        f"{viscosity:.2e} Pa·s, and its density (ρ) is {density:.3f} kg/m³.\n\n"
        f"Calculate the kinematic viscosity (ν) of {fluid_name} in units of **{target_units}**."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Dynamic Viscosity (μ): {viscosity:.2e} Pa·s\n"
        f"- Density (ρ): {density:.3f} kg/m³\n\n"
        
        f"**Step 1:** State the formula for kinematic viscosity.\n"
        f"Kinematic viscosity (ν) is defined as the ratio of dynamic viscosity to density:\n"
        f"ν = μ / ρ\n\n"
        f"\n\n"
        
        f"**Step 2:** Calculate the kinematic viscosity in SI units (m²/s).\n"
        f"Substitute the given values into the formula:\n"
        f"ν = ({viscosity:.2e} Pa·s) / ({density:.3f} kg/m³)\n"
        f"ν = {nu_si:.3e} m²/s\n\n"
    )
    
    # Add the unit conversion step only if needed
    if "Stokes" in target_units:
        # Fixed formatting here to avoid 0.0000 for small values
        solution += (
            f"**Step 3:** Convert the result to Stokes (St).\n"
            f"The conversion factor is 1 m²/s = 10,000 St.\n"
            f"ν = ({nu_si:.3e} m²/s) * (10,000 St / 1 m²/s)\n"
            f"ν = {final_nu:.4e} St\n\n"
        )
        
    solution += (
        f"**Answer:**\n"
        # Fixed formatting here to avoid 0.0000 for small values
        f"The kinematic viscosity of {fluid_name} is **{final_nu:.4e} {target_units}**."
    )

    return question, solution


# Template 3 (Intermediate)
def template_gas_viscosity_kinetic_theory():
    """
    Gas Viscosity Estimation from Kinetic Theory

    Scenario:
        This template uses the Chapman-Enskog equation, derived from the kinetic
        theory of gases, to estimate the dynamic viscosity (μ) of a low-density gas.
        This model connects a macroscopic property (viscosity) to molecular
        properties (molar mass, size).

        The relevant equation is:
            μ = (2.6693e-6 * sqrt(M*T)) / (σ² * Ω_μ)
        
        Where:
            - μ = viscosity in Pa·s
            - M = Molar Mass in g/mol
            - T = Absolute Temperature in K
            - σ = Lennard-Jones molecular diameter in Angstroms (Å)
            - Ω_μ = Collision integral (dimensionless)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to estimate the gas viscosity.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    gas_name, (molar_mass, sigma, epsilon) = random.choice(list(GAS_MOLECULAR_PARAMS.items()))
    
    # Temperature in Kelvin
    temperature_K = random.randint(250, 600)
    
    # Collision integral (dimensionless), kept close to 1.0 for simplicity
    omega_mu = round(random.uniform(0.95, 1.05), 3)

    # 2. Perform the core calculation
    
    # Numerator of the Chapman-Enskog equation
    numerator = 2.6693e-6 * math.sqrt(molar_mass * temperature_K)
    
    # Denominator of the Chapman-Enskog equation
    denominator = (sigma**2) * omega_mu
    
    # Final viscosity calculation
    viscosity = numerator / denominator

    # 3. Generate the question and solution strings
    
    question = (
        f"Estimate the dynamic viscosity (μ) of {gas_name} gas at a temperature of {temperature_K} K.\n\n"
        f"Use the Chapman-Enskog equation for low-density gases:\n"
        f"μ = (2.6693e-6 * sqrt(M*T)) / (σ^2 * Ω_μ)\n\n"
        f"You are given the following parameters for {gas_name}:\n"
        f"- Molar Mass (M) = {molar_mass} g/mol\n"
        f"- Lennard-Jones diameter (σ) = {sigma} Å\n"
        f"- Collision Integral (Ω_μ) = {omega_mu}\n\n"
        f"Provide your answer in units of Pascal-seconds (Pa·s)."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Gas: {gas_name}\n"
        f"- Temperature (T): {temperature_K} K\n"
        f"- Molar Mass (M): {molar_mass} g/mol\n"
        f"- Lennard-Jones diameter (σ): {sigma} Å\n"
        f"- Collision Integral (Ω_μ): {omega_mu}\n\n"
        
        f"**Step 1:** State the Chapman-Enskog equation.\n"
        f"The formula for estimating the viscosity of a low-density gas is:\n"
        f"μ = (2.6693e-6 * sqrt(M*T)) / (σ^2 * Ω_μ)\n\n"
        
        f"**Step 2:** Substitute the given values into the equation.\n"
        f"μ = (2.6693e-6 * sqrt({molar_mass} * {temperature_K})) / (({sigma})^2 * {omega_mu})\n\n"
        
        f"**Step 3:** Calculate the numerator and denominator.\n"
        f"- Numerator = 2.6693e-6 * sqrt({molar_mass * temperature_K:.2f}) = {numerator:.4e}\n"
        f"- Denominator = {sigma**2:.3f} * {omega_mu} = {denominator:.3f}\n\n"
        
        f"**Step 4:** Calculate the final viscosity.\n"
        f"μ = {numerator:.4e} / {denominator:.3f}\n"
        f"μ = {viscosity:.3e} Pa·s\n\n"
        
        f"**Answer:**\n"
        f"The estimated dynamic viscosity of {gas_name} at {temperature_K} K is **{viscosity:.3e} Pa·s**."
    )

    return question, solution


# Template 4 (Intermediate)
def template_reynolds_number_flow_regime():
    """
    Reynolds Number and Flow Regime Calculation

    Scenario:
        This template introduces the Reynolds number (Re), the dimensionless ratio of
        inertial forces to viscous forces in a fluid. Its value is critical for
        predicting whether a flow is smooth (laminar) or chaotic (turbulent).
        The calculation depends on a characteristic length, which varies with geometry.

        The relevant equation is:
            Re = (ρ * v * L) / μ
        
        Where:
            - ρ = density
            - v = average velocity
            - L = characteristic length (e.g., pipe diameter or plate length)
            - μ = dynamic viscosity

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the Reynolds number and find the flow regime.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    all_fluids = {**COMMON_LIQUIDS, **COMMON_GASES}
    fluid_name, (density, viscosity) = random.choice(list(all_fluids.items()))
    
    # Randomly choose a flow geometry
    geometry = random.choice(['pipe', 'flat plate'])
    
    if geometry == 'pipe':
        # Pipe flow parameters
        characteristic_length = round(random.uniform(0.01, 0.5), 3) # Diameter in m
        velocity = round(random.uniform(0.1, 5.0), 2)
        length_symbol = "D"
        length_name = "diameter"
        scenario_description = f"{fluid_name} is flowing through a smooth circular pipe with an internal {length_name} of {characteristic_length} m."
        
    else: # flat plate
        # Flat plate flow parameters
        characteristic_length = round(random.uniform(0.1, 2.0), 2) # Length in m
        velocity = round(random.uniform(1.0, 20.0), 1)
        length_symbol = "L"
        length_name = "length"
        scenario_description = f"A flow of {fluid_name} moves over a smooth, thin flat plate with a {length_name} of {characteristic_length} m."

    # 2. Perform the core calculation
    reynolds_number = (density * velocity * characteristic_length) / viscosity

    # Determine flow regime based on geometry
    if geometry == 'pipe':
        if reynolds_number < 2300:
            regime = "laminar"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is less than 2300, the flow is in the **laminar** regime."
        elif reynolds_number < 4000:
            regime = "transitional"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is between 2300 and 4000, the flow is in the **transitional** regime."
        else:
            regime = "turbulent"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is greater than 4000, the flow is in the **turbulent** regime."
    else: # flat plate
        if reynolds_number < 500000:
            regime = "laminar"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is less than the critical value of 500,000, the flow is considered **laminar**."
        else:
            regime = "turbulent"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is greater than the critical value of 500,000, the flow is considered **turbulent**."

    # 3. Generate the question and solution strings
    question = (
        f"{scenario_description} The average velocity of the flow is {velocity} m/s.\n\n"
        f"The properties of {fluid_name} are:\n"
        f"- Density (ρ) = {density} kg/m³\n"
        f"- Dynamic Viscosity (μ) = {viscosity:.2e} Pa·s\n\n"
        f"Based on this information:\n"
        f"a) Calculate the Reynolds number (Re).\n"
        f"b) Determine the flow regime (laminar, transitional, or turbulent)."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Density (ρ): {density} kg/m³\n"
        f"- Dynamic Viscosity (μ): {viscosity:.2e} Pa·s\n"
        f"- Average Velocity (v): {velocity} m/s\n"
        f"- Characteristic Length ({length_symbol}): {characteristic_length} m ({length_name} of the {geometry})\n\n"
        
        f"**Step 1:** State the formula for the Reynolds number.\n"
        f"The Reynolds number is the ratio of inertial forces to viscous forces.\n"
        f"Re = (ρ * v * {length_symbol}) / μ\n\n"
        
        f"**Step 2:** Substitute the values and calculate Re.\n"
        f"Re = ({density} * {velocity} * {characteristic_length}) / {viscosity:.2e}\n"
        f"Re = {reynolds_number:,.0f}\n\n"
        
        f"**Step 3:** Determine the flow regime.\n"
        f"For flow in a **{geometry}**, we compare the calculated Re to the standard critical values.\n"
        f"{regime_explanation}\n\n"

        f"**Answer:**\n"
        f"a) The Reynolds number is approximately **{reynolds_number:,.0f}**.\n"
        f"b) The flow regime is **{regime}**."
    )

    return question, solution


# Template 5 (Advanced)
def template_power_law_fluid_shear():
    """
    Shear Stress Calculation for a Power-Law Fluid

    Scenario:
        This template extends Newton's Law of Viscosity to non-Newtonian fluids
        that follow the Ostwald-de Waele (Power Law) model. Unlike Newtonian fluids,
        their viscosity is dependent on the shear rate. This model is crucial for
        describing common substances like paint, ketchup, and suspensions.

        The relevant equations are:
            - Shear Stress: τ_yx = K * (dvx/dy)^n
            - Apparent Viscosity: η = K * |dvx/dy|^(n-1)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for shear stress and apparent viscosity.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    fluid_name, (K, n) = random.choice(list(POWER_LAW_FLUIDS.items()))
    
    # Velocity of the moving plate in m/s
    V = round(random.uniform(0.1, 2.0), 2)
    
    # Distance between plates in cm, converted to meters for calculation
    Y_cm = round(random.uniform(0.5, 5.0), 2)
    Y_m = Y_cm / 100

    # 2. Perform the core calculation
    velocity_gradient = abs(V / Y_m)
    shear_stress = K * (velocity_gradient ** n)
    apparent_viscosity = K * (velocity_gradient ** (n - 1))
    
    # Determine the fluid behavior for the explanation
    if n < 1:
        behavior = f"shear-thinning (pseudoplastic), because its power-law index n ({n}) is less than 1."
    elif n > 1:
        behavior = f"shear-thickening (dilatant), because its power-law index n ({n}) is greater than 1."
    else:
        behavior = "Newtonian, because its power-law index n is exactly 1."

    # 3. Generate the question and solution strings
    question = (
        f"A non-Newtonian fluid, {fluid_name.lower()}, is placed between two parallel plates "
        f"separated by {Y_cm} cm. The top plate moves at a constant velocity of {V} m/s, "
        f"creating a linear velocity profile in the fluid.\n\n"
        f"The fluid follows the power-law model with the following parameters:\n"
        f"- Consistency Index (K) = {K} Pa·s^n\n"
        f"- Power-Law Index (n) = {n}\n\n"
        f"Calculate:\n"
        f"a) The shear stress (τ_yx) on the fluid.\n"
        f"b) The apparent viscosity (η) at this shear rate."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Consistency Index (K): {K} Pa·s^n\n"
        f"- Power-Law Index (n): {n}\n"
        f"- Plate Velocity (V): {V} m/s\n"
        f"- Plate Separation (Y): {Y_cm} cm = {Y_m} m\n\n"
        
        f"**Step 1:** Characterize the Fluid Behavior.\n"
        f"The fluid is {behavior} This means its apparent viscosity will change with the rate of shear.\n\n"
        
        f"**Step 2:** Calculate the Velocity Gradient (Shear Rate).\n"
        f"Assuming a linear velocity profile:\n"
        f"Shear Rate (dvx/dy) = V / Y = {V} m/s / {Y_m} m = {velocity_gradient:.2f} 1/s\n\n"
        
        f"**Step 3:** Calculate the Shear Stress (τ_yx).\n"
        f"Using the power-law formula: τ_yx = K * (dvx/dy)^n\n"
        f"τ_yx = {K} * ({velocity_gradient:.2f})^{n}\n"
        f"τ_yx = {shear_stress:.3f} Pa\n\n"
        
        f"**Step 4:** Calculate the Apparent Viscosity (η).\n"
        f"Apparent viscosity is the effective viscosity at a specific shear rate.\n"
        f"η = K * |dvx/dy|^(n-1)\n"
        f"η = {K} * ({velocity_gradient:.2f})^({n - 1})\n"
        f"η = {apparent_viscosity:.4f} Pa·s\n\n"
        
        f"**Answer:**\n"
        f"a) The shear stress on the fluid is **{shear_stress:.3f} Pa**.\n"
        f"b) The apparent viscosity at this shear rate is **{apparent_viscosity:.4f} Pa·s**."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each viscosity and mechanisms of momentum transport template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/transport_phenomena/viscosity_and_momentum_transport.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_newtons_law_shear_stress, "newtons_law_shear_stress", "Easy"),
        (template_kinematic_viscosity, "kinematic_viscosity", "Easy"),
        (template_gas_viscosity_kinetic_theory, "gas_viscosity_kinetic_theory", "Intermediate"),
        (template_reynolds_number_flow_regime, "reynolds_number_flow_regime", "Intermediate"),
        (template_power_law_fluid_shear, "power_law_fluid_shear", "Advanced"),
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
                "domain": "transport_phenomena",
                "area": "viscosity_and_momentum_transport",
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
