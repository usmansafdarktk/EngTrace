import random
import math
from data.templates.branches.chemical_engineering.constants import COMMON_LIQUIDS, GRAVITATIONAL_ACCELERATION


# Template 1 (Easy)
def template_falling_film_max_velocity():
    """
    Falling Film Maximum Velocity Calculation

    Scenario:
        This template calculates the maximum velocity of a liquid film flowing
        down an inclined plane under gravity. The maximum velocity occurs at the
        free surface (the liquid-air interface), where shear stress is zero.
        The calculation is based on the final velocity profile equation derived
        from a shell momentum balance for this system.

        Core Equation:
            v_z_max = (rho * g * delta^2 * cos(beta)) / (2 * mu)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the maximum film velocity.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    # Choose a random fluid name (key) from the dictionary
    fluid_name = random.choice(list(COMMON_LIQUIDS.keys()))
    # Get the corresponding properties (the tuple of density and viscosity)
    fluid_density, fluid_viscosity = COMMON_LIQUIDS[fluid_name]
    
    # Film thickness in mm
    film_thickness_mm = round(random.uniform(0.5, 5.0), 2)
    
    # Angle of inclination from the vertical in degrees
    inclination_angle_deg = random.randint(5, 85)
    
    # Constants
    g = GRAVITATIONAL_ACCELERATION

    # 2. Perform unit conversions and core calculation
    # Convert film thickness from mm to m
    film_thickness_m = film_thickness_mm / 1000.0
    
    # Convert inclination angle from degrees to radians for math functions
    inclination_angle_rad = math.radians(inclination_angle_deg)
    
    # Calculate the maximum velocity
    v_max = (fluid_density * g * (film_thickness_m**2) * math.cos(inclination_angle_rad)) / (2 * fluid_viscosity)

    # 3. Generate the question and solution strings
    question = (
        f"A thin film of {fluid_name} is flowing down a flat surface inclined at an "
        f"angle of {inclination_angle_deg} degrees from the vertical. The film thickness is "
        f"measured to be {film_thickness_mm} mm.\n\n"
        f"Given the fluid's density is {fluid_density} kg/m^3 and its viscosity is "
        f"{fluid_viscosity} Pa·s, calculate the maximum velocity of the film (at the "
        f"liquid-air interface). Assume the flow is laminar.\n\n"
        f"Use a value of g = {g} m/s^2."
    )

    solution = (
        f"**Step 1:** Identify Given Information\n"
        f"First, we list the known values from the problem statement.\n"
        f"- Fluid: {fluid_name}\n"
        f"- Fluid Density (rho): {fluid_density} kg/m^3\n"
        f"- Fluid Viscosity (mu): {fluid_viscosity} Pa·s\n"
        f"- Film Thickness (delta): {film_thickness_mm} mm\n"
        f"- Inclination Angle (beta): {inclination_angle_deg} degrees\n"
        f"- Gravitational Acceleration (g): {g} m/s^2\n\n"
        
        f"**Step 2:** Perform Unit Conversions\n"
        f"The core equation requires all units to be in the SI base system. We must convert the film thickness from millimeters to meters.\n"
        f"delta = {film_thickness_mm} mm * (1 m / 1000 mm) = {film_thickness_m} m\n\n"
        
        f"**Step 3:** State the Core Equation\n"
        f"For a laminar falling film, the maximum velocity (v_z_max) at the liquid-air interface is given by the equation derived from the shell momentum balance:\n"
        f"v_z_max = (rho * g * delta^2 * cos(beta)) / (2 * mu)\n\n"
        
        f"**Step 4:** Substitute Values into the Equation\n"
        f"Now, we substitute our known values (with correct units) into the equation.\n"
        f"v_z_max = (({fluid_density} kg/m^3) * ({g} m/s^2) * ({film_thickness_m} m)^2 * cos({inclination_angle_deg} deg)) / (2 * ({fluid_viscosity} Pa·s))\n\n"
        
        f"**Step 5:** Calculate the Final Velocity\n"
        f"Performing the calculation gives the maximum velocity.\n"
        f"v_z_max = {round(v_max, 5)} m/s\n\n"
        
        f"**Answer:** The maximum velocity of the {fluid_name} film is {round(v_max, 5)} m/s."
    )

    return question, solution


# Template 2 (Easy)
def template_hagen_poiseuille_flowrate():
    """
    Hagen-Poiseuille Volumetric Flow Rate Calculation

    Scenario:
        This template calculates the volumetric flow rate (Q) for a fluid
        experiencing laminar flow through a straight, circular pipe of
        constant cross-section. It's based on the exact solution for
        pressure-driven flow derived from a shell momentum balance.

        Core Equation (Hagen-Poiseuille Equation):
            Q = (pi * (P0 - PL) * R^4) / (8 * mu * L)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the volumetric flow rate.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    # Choose a random fluid name (key) from the dictionary
    fluid_name = random.choice(list(COMMON_LIQUIDS.keys()))
    # Get the corresponding properties (the tuple of density and viscosity)
    fluid_density, fluid_viscosity_Pas = COMMON_LIQUIDS[fluid_name]

    # Pipe radius in cm
    pipe_radius_cm = round(random.uniform(0.5, 5.0), 2)
    # Pipe length in m
    pipe_length_m = round(random.uniform(5.0, 100.0), 1)
    # Pressure drop in kPa
    pressure_drop_kPa = random.randint(50, 500)
    # Convert viscosity to cP for the problem statement (1 Pa·s = 1000 cP)
    fluid_viscosity_cP = fluid_viscosity_Pas * 1000

    # 2. Perform unit conversions and core calculation
    # Convert radius from cm to m
    pipe_radius_m = pipe_radius_cm / 100.0
    # Convert pressure drop from kPa to Pa
    pressure_drop_Pa = pressure_drop_kPa * 1000

    # Calculate the volumetric flow rate (Q)
    # Q = (pi * delta_P * R^4) / (8 * mu * L)
    q_flow_rate = (math.pi * pressure_drop_Pa * (pipe_radius_m**4)) / (8 * fluid_viscosity_Pas * pipe_length_m)

    # 3. Generate the question and solution strings
    question = (
        f"A fluid, {fluid_name}, is flowing through a smooth circular pipe. "
        f"The pipe has an inner radius of {pipe_radius_cm} cm and a total length of {pipe_length_m} m. "
        f"The pressure drop across the length of the pipe is measured to be {pressure_drop_kPa} kPa.\n\n"
        f"Given that the viscosity of {fluid_name} is approximately {fluid_viscosity_cP:.2f} cP, "
        f"calculate the volumetric flow rate (Q) in m^3/s. Assume the flow is laminar."
    )

    solution = (
        f"**Step 1:** Identify Given Information\n"
        f"First, we list the known values from the problem statement.\n"
        f"- Fluid: {fluid_name}\n"
        f"- Pipe Radius (R): {pipe_radius_cm} cm\n"
        f"- Pipe Length (L): {pipe_length_m} m\n"
        f"- Pressure Drop (P0 - PL): {pressure_drop_kPa} kPa\n"
        f"- Fluid Viscosity (mu): {fluid_viscosity_cP:.2f} cP\n\n"

        f"**Step 2:** Perform Unit Conversions\n"
        f"The Hagen-Poiseuille equation requires all units to be in the SI base system. We must convert the radius, pressure drop, and viscosity.\n"
        f"- Radius: R = {pipe_radius_cm} cm * (1 m / 100 cm) = {pipe_radius_m} m\n"
        f"- Pressure Drop: P0 - PL = {pressure_drop_kPa} kPa * (1000 Pa / 1 kPa) = {pressure_drop_Pa} Pa\n"
        f"- Viscosity: mu = {fluid_viscosity_cP:.2f} cP * (1 Pa·s / 1000 cP) = {fluid_viscosity_Pas} Pa·s\n\n"

        f"**Step 3:** State the Core Equation\n"
        f"The volumetric flow rate (Q) for laminar flow in a circular pipe is given by the Hagen-Poiseuille equation:\n"
        f"Q = (pi * (P0 - PL) * R^4) / (8 * mu * L)\n\n"

        f"**Step 4:** Substitute Values into the Equation\n"
        f"Now, we substitute the converted SI values into the equation.\n"
        f"Q = (pi * ({pressure_drop_Pa} Pa) * ({pipe_radius_m} m)^4) / (8 * ({fluid_viscosity_Pas} Pa·s) * ({pipe_length_m} m))\n\n"

        f"**Step 5:** Calculate the Final Flow Rate\n"
        f"Performing the calculation gives the volumetric flow rate.\n"
        f"Q = {q_flow_rate:.6f} m^3/s\n\n"

        f"**Answer:** The volumetric flow rate of {fluid_name} through the pipe is {q_flow_rate:.6f} m^3/s."
    )

    return question, solution


# Template 3 (Intermediate)
def template_annulus_flowrate():
    """
    Volumetric Flow Rate in an Annulus Calculation

    Scenario:
        This template calculates the volumetric flow rate (Q) for a fluid in
        laminar flow through the concentric cylindrical region between two pipes (an annulus).
        The equation is derived from the shell momentum balance for this specific geometry.

        Core Equation:
            Q = (pi*(P0-PL)*R^4)/(8*mu*L) * [(1-kappa^4) - ((1-kappa^2)^2 / ln(1/kappa))]

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the volumetric flow rate in an annulus.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    # Choose a random fluid name (key) from the dictionary
    fluid_name = random.choice(list(COMMON_LIQUIDS.keys()))
    # Get the corresponding properties (the tuple of density and viscosity)
    fluid_density, fluid_viscosity_Pas = COMMON_LIQUIDS[fluid_name]

    # Define radii and ensure inner radius is smaller than outer radius
    outer_radius_cm = round(random.uniform(2.0, 10.0), 2)
    kappa = round(random.uniform(0.2, 0.8), 2) # Ratio of inner to outer radius
    inner_radius_cm = round(kappa * outer_radius_cm, 2)

    # Pipe length in m
    pipe_length_m = round(random.uniform(10.0, 150.0), 1)
    # Pressure drop in kPa
    pressure_drop_kPa = random.randint(100, 1000)

    # 2. Perform unit conversions and core calculation
    # Convert radii from cm to m
    outer_radius_m = outer_radius_cm / 100.0
    inner_radius_m = inner_radius_cm / 100.0
    
    # Convert pressure drop from kPa to Pa
    pressure_drop_Pa = pressure_drop_kPa * 1000

    # The dimensionless ratio, kappa
    kappa_val = inner_radius_m / outer_radius_m

    # Calculate the flow rate using the annulus equation
    term1 = (math.pi * pressure_drop_Pa * (outer_radius_m**4)) / (8 * fluid_viscosity_Pas * pipe_length_m)
    shape_factor = (1 - kappa_val**4) - (((1 - kappa_val**2)**2) / math.log(1 / kappa_val))
    q_flow_rate = term1 * shape_factor

    # 3. Generate the question and solution strings
    question = (
        f"Laminar flow of {fluid_name} occurs in the annular space between two concentric pipes. "
        f"The inner pipe has an outer radius of {inner_radius_cm} cm, and the outer pipe has an inner radius of {outer_radius_cm} cm. "
        f"The concentric pipes have a length of {pipe_length_m} m.\n\n"
        f"A pressure drop of {pressure_drop_kPa} kPa is maintained over the length of the pipes. "
        f"The fluid viscosity is {fluid_viscosity_Pas} Pa·s.\n\n"
        f"Calculate the volumetric flow rate (Q) through the annular space in m^3/s."
    )

    solution = (
        f"**Step 1:** Identify Given Information\n"
        f"First, we list the known values from the problem statement.\n"
        f"- Fluid: {fluid_name}\n"
        f"- Inner Radius (R_inner): {inner_radius_cm} cm\n"
        f"- Outer Radius (R_outer): {outer_radius_cm} cm\n"
        f"- Pipe Length (L): {pipe_length_m} m\n"
        f"- Pressure Drop (P0 - PL): {pressure_drop_kPa} kPa\n"
        f"- Fluid Viscosity (mu): {fluid_viscosity_Pas} Pa·s\n\n"

        f"**Step 2:** Perform Unit Conversions\n"
        f"The equation requires all units to be in the SI base system. We must convert the radii and the pressure drop.\n"
        f"- Inner Radius: R_inner = {inner_radius_cm} cm * (1 m / 100 cm) = {inner_radius_m} m\n"
        f"- Outer Radius: R_outer = {outer_radius_cm} cm * (1 m / 100 cm) = {outer_radius_m} m\n"
        f"- Pressure Drop: P0 - PL = {pressure_drop_kPa} kPa * (1000 Pa / 1 kPa) = {pressure_drop_Pa} Pa\n\n"
        
        f"**Step 3:** Calculate the Dimensionless Ratio (kappa)\n"
        f"Kappa (k) is the ratio of the inner radius to the outer radius.\n"
        f"kappa = R_inner / R_outer = {inner_radius_m} m / {outer_radius_m} m = {kappa_val:.3f}\n\n"

        f"**Step 4:** State the Core Equation\n"
        f"The volumetric flow rate (Q) for laminar flow in an annulus is:\n"
        f"Q = (pi * (P0 - PL) * R_outer^4) / (8 * mu * L) * [ (1 - kappa^4) - ((1 - kappa^2)^2 / ln(1/kappa)) ]\n\n"

        f"**Step 5:** Substitute Values and Calculate\n"
        f"We substitute the converted SI values into the equation. Let's calculate the main term and the shape factor separately for clarity.\n"
        f"Main Term = (pi * {pressure_drop_Pa} Pa * ({outer_radius_m} m)^4) / (8 * {fluid_viscosity_Pas} Pa·s * {pipe_length_m} m) = {term1:.6f}\n"
        f"Shape Factor = [ (1 - {kappa_val:.3f}^4) - ((1 - {kappa_val:.3f}^2)^2 / ln(1/{kappa_val:.3f})) ] = {shape_factor:.4f}\n\n"
        f"Q = Main Term * Shape Factor\n"
        f"Q = {term1:.6f} * {shape_factor:.4f} = {q_flow_rate:.7f} m^3/s\n\n"

        f"Answer: The volumetric flow rate of {fluid_name} through the annular space is {q_flow_rate:.7f} m^3/s."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each shell momentum balances and velocity distributions in 
    laminar flow template with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/transport_phenomena/shell_momentum_balances.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_falling_film_max_velocity, "falling_film_max_velocity", "Easy"),
        (template_hagen_poiseuille_flowrate, "hagen_poiseuille_flowrate", "Easy"),
        (template_annulus_flowrate, "annulus_flowrate", "Intermediate"),
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
