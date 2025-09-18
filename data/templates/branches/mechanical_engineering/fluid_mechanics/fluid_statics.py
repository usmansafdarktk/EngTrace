import random
import math
from data.templates.branches.mechanical_engineering.constants import GRAVITY, ATMOSPHERIC_PRESSURE_KPA, FLUID_DENSITIES, MATERIAL_DENSITIES, OBJECT_SHAPES, OBJECT_MATERIALS, PIPE_FLUIDS, MANOMETER_FLUIDS


# Template 1 (Easy)
def template_hydrostatic_pressure_at_depth():
    """
    Fluid Statics: Hydrostatic Pressure at Depth

    Scenario:
        This template generates a fundamental problem to calculate the absolute
        pressure at a specified depth within a fluid. It considers both the
        pressure exerted by the fluid column and the pressure at the free surface,
        applying the basic hydrostatic pressure equation.

    Core Equation:
        p_absolute = p_surface + (rho * g * h)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the absolute pressure at a certain depth.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values

    # Randomly select a fluid and its properties
    fluid_name, density_rho = random.choice(list(FLUID_DENSITIES.items()))

    # Randomize depth in meters
    depth_h = round(random.uniform(5.0, 500.0), 2)

    # Randomly decide if the surface pressure is atmospheric or a specified gauge pressure
    is_surface_atmospheric = random.choice([True, False])

    if is_surface_atmospheric:
        surface_pressure_kpa = ATMOSPHERIC_PRESSURE_KPA
        pressure_type_desc = "is atmospheric pressure"
    else:
        # Generate a random gauge pressure at the surface
        surface_pressure_kpa = round(random.uniform(10.0, 200.0), 2)
        pressure_type_desc = f"has a gauge pressure of {surface_pressure_kpa} kPa"

    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculations for the solution

    # Step A: Ensure all units are consistent (convert surface pressure to Pascals)
    surface_pressure_pa = surface_pressure_kpa * 1000

    # Step B: Calculate the pressure increase due to the fluid column (Gauge Pressure)
    # This is the result of p_gauge = rho * g * h
    pressure_increase_pa = density_rho * GRAVITY * depth_h

    # Step C: Calculate the final absolute pressure in Pascals
    absolute_pressure_pa = surface_pressure_pa + pressure_increase_pa

    # Step D: Convert the final answer back to kilopascals (kPa) for readability
    absolute_pressure_kpa = absolute_pressure_pa / 1000

    # 3. Generate the question and solution strings

    question = (
        f"An object is located {depth_h} m below the surface of a tank containing "
        f"{fluid_name.lower()}. The pressure at the free surface {pressure_type_desc}. "
        f"Assuming the density of {fluid_name.lower()} is {density_rho} kg/m^3, "
        f"what is the absolute pressure at this depth?"
    )

    solution = (
        f"**Given:**\n"
        f"Fluid: {fluid_name}\n"
        f"Density of Fluid (rho): {density_rho} kg/m^3\n"
        f"Depth (h): {depth_h} m\n"
        f"Surface Pressure (p_surface): {surface_pressure_kpa} kPa\n"
        f"Acceleration due to Gravity (g): {GRAVITY} m/s^2\n\n"

        f"**Step 1:** Ensure Consistent Units\n"
        f"The calculations require pressure to be in Pascals (Pa) to be consistent with other SI units.\n"
        f"p_surface = {surface_pressure_kpa} kPa * 1000 = {surface_pressure_pa:.2f} Pa\n\n"

        f"**Step 2:** Calculate the Pressure Increase Due to the Fluid Column\n"
        f"The pressure exerted by the fluid at a given depth is calculated using the formula: p_increase = rho * g * h.\n"
        f"p_increase = {density_rho} kg/m^3 * {GRAVITY} m/s^2 * {depth_h} m\n"
        f"p_increase = {pressure_increase_pa:.2f} Pa\n\n"

        f"**Step 3:** Calculate the Absolute Pressure at the Specified Depth\n"
        f"The absolute pressure is the sum of the surface pressure and the pressure increase due to the fluid column.\n"
        f"Formula: p_absolute = p_surface + p_increase\n"
        f"p_absolute = {surface_pressure_pa:.2f} Pa + {pressure_increase_pa:.2f} Pa\n"
        f"p_absolute = {absolute_pressure_pa:.2f} Pa\n\n"

        f"**Step 4:** Convert the Final Answer to Kilopascals (kPa)\n"
        f"For convenience, we convert the final pressure from Pascals back to kilopascals.\n"
        f"p_absolute = {absolute_pressure_pa:.2f} Pa / 1000 = {round(absolute_pressure_kpa, precision)} kPa\n\n"

        f"**Answer:**\n"
        f"The absolute pressure at a depth of {depth_h} m is {round(absolute_pressure_kpa, precision)} kPa."
    )

    return question, solution


# Template 2 (Easy)
def template_basic_buoyant_force():
    """
    Fluid Statics: Basic Buoyant Force on a Submerged Object

    Scenario:
        This template applies Archimedes' principle to calculate the buoyant
        force on a fully submerged object with a known volume. It reinforces
        the direct relationship between buoyant force and the weight of the
        displaced fluid.

    Core Equation:
        F_buoyant = rho_fluid * g * V_displaced

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the buoyant force on an object.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values

    # Randomly select a fluid and its properties
    fluid_name, density_rho = random.choice(list(FLUID_DENSITIES.items()))

    # Randomly select descriptive properties for the object
    shape = random.choice(OBJECT_SHAPES)
    material = random.choice(OBJECT_MATERIALS)

    # Randomize the object's volume in cubic meters
    object_volume = round(random.uniform(0.05, 2.5), 3)

    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculations for the solution

    # Step A: Calculate the buoyant force in Newtons (N)
    # Since the object is fully submerged, the displaced volume equals the object's volume.
    buoyant_force_n = density_rho * GRAVITY * object_volume

    # Step B: Convert to kilonewtons (kN) if the value is large, for better readability
    buoyant_force_kn = buoyant_force_n / 1000

    # 3. Generate the question and solution strings

    question = (
        f"A solid {material} {shape} with a total volume of {object_volume} m^3 is "
        f"fully submerged in a tank of {fluid_name.lower()}. "
        f"Given that the density of {fluid_name.lower()} is {density_rho} kg/m^3, "
        f"calculate the buoyant force acting on the {shape}."
    )

    solution = (
        f"**Given:**\n"
        f"Object Volume (V): {object_volume} m^3\n"
        f"Fluid: {fluid_name}\n"
        f"Density of Fluid (rho): {density_rho} kg/m^3\n"
        f"Acceleration due to Gravity (g): {GRAVITY} m/s^2\n\n"

        f"**Step 1:** Identify the Principle\n"
        f"According to Archimedes' principle, the buoyant force (F_buoyant) on a submerged "
        f"object is equal to the weight of the fluid it displaces.\n"
        f"Since the object is fully submerged, the volume of displaced fluid is equal to the "
        f"volume of the object itself.\n\n"

        f"**Step 2:** Apply the Buoyant Force Formula\n"
        f"The formula for buoyant force is: F_buoyant = rho_fluid * g * V_displaced\n"
        f"F_buoyant = {density_rho} kg/m^3 * {GRAVITY} m/s^2 * {object_volume} m^3\n"
        f"F_buoyant = {round(buoyant_force_n, precision)} N\n\n"

        f"**Step 3:** Convert the Result to Kilonewtons (kN) (Optional)\n"
        f"For larger force values, it is common to express the result in kilonewtons.\n"
        f"F_buoyant = {round(buoyant_force_n, precision)} N / 1000 = {round(buoyant_force_kn, precision)} kN\n\n"

        f"**Answer:**\n"
        f"The buoyant force acting on the {shape} is {round(buoyant_force_n, precision)} N, "
        f"which is equivalent to {round(buoyant_force_kn, precision)} kN."
    )

    return question, solution


# Template 3 (Intermediate)
def template_utube_manometer():
    """
    Fluid Statics: U-Tube Manometer Pressure Measurement

    Scenario:
        This template generates a classic problem involving a U-tube manometer
        used to measure the gauge pressure of a fluid inside a pipe. The solution
        requires applying the principle of hydrostatic equilibrium by balancing
        the pressure contributions from different fluid columns.

    Core Equation:
        P_gauge = (rho_2 * g * h2) - (rho_1 * g * h1)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the gauge pressure in a pipe.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values

    # Randomly select a fluid for the pipe
    pipe_fluid_name, rho1 = random.choice(list(PIPE_FLUIDS.items()))

    # Randomly select a fluid for the manometer
    manometer_fluid_name, rho2 = random.choice(list(MANOMETER_FLUIDS.items()))

    # --- Ensure physical realism: manometer fluid must be denser ---
    # If the chosen pipe fluid is denser than the manometer fluid, re-select
    # the manometer fluid until it is denser. This avoids nonsensical scenarios.
    while rho2 <= rho1:
        manometer_fluid_name, rho2 = random.choice(list(MANOMETER_FLUIDS.items()))

    # Randomize the vertical heights, in meters
    # h1: distance from pipe centerline down to the fluid interface
    h1 = round(random.uniform(0.1, 0.5), 2)
    # h2: the height difference between the two manometer fluid columns
    h2 = round(random.uniform(0.05, 0.75), 2)

    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculations for the solution

    # Step A: Calculate the pressure contribution from the pipe fluid column (rho1*g*h1)
    pressure_term1_pa = rho1 * GRAVITY * h1

    # Step B: Calculate the pressure contribution from the manometer fluid column (rho2*g*h2)
    pressure_term2_pa = rho2 * GRAVITY * h2

    # Step C: Calculate the gauge pressure in Pascals (Pa)
    gauge_pressure_pa = pressure_term2_pa - pressure_term1_pa

    # Step D: Convert the final answer to kilopascals (kPa) for readability
    gauge_pressure_kpa = gauge_pressure_pa / 1000

    # 3. Generate the question and solution strings

    question = (
        f"A U-tube manometer using {manometer_fluid_name.lower()} (density = {rho2} kg/m^3) is "
        f"connected to a pipe carrying {pipe_fluid_name.lower()} (density = {rho1} kg/m^3). "
        f"The interface between the two fluids is {h1} m below the centerline of the pipe. "
        f"The level of {manometer_fluid_name.lower()} in the arm open to the atmosphere is {h2} m "
        f"higher than the interface. What is the gauge pressure in the pipe?"
    )

    solution = (
        f"**Given:**\n"
        f"Pipe Fluid: {pipe_fluid_name}, Density (rho1) = {rho1} kg/m^3\n"
        f"Manometer Fluid: {manometer_fluid_name}, Density (rho2) = {rho2} kg/m^3\n"
        f"Height from pipe centerline to interface (h1): {h1} m\n"
        f"Height difference in manometer fluid (h2): {h2} m\n"
        f"Acceleration due to Gravity (g): {GRAVITY} m/s^2\n\n"

        f"**Step 1:** State the Principle of Manometry\n"
        f"We can determine the pressure in the pipe by starting at the pipe's centerline, moving "
        f"through the fluid columns to the open end, and balancing the pressures. The pressure at the "
        f"same level within a continuous fluid at rest is equal.\n\n"

        f"**Step 2:** Formulate the Pressure Balance Equation\n"
        f"Let's establish a pressure equation starting from the pipe (P_pipe) and ending at the atmosphere (P_atm). "
        f"The reference level is the interface between the pipe fluid and the manometer fluid.\n"
        f"Pressure from pipe side at interface: P_pipe + (rho1 * g * h1)\n"
        f"Pressure from atmosphere side at interface: P_atm + (rho2 * g * h2)\n"
        f"Equating these gives: P_pipe + (rho1 * g * h1) = P_atm + (rho2 * g * h2)\n\n"

        f"**Step 3:** Solve for the Gauge Pressure\n"
        f"Gauge pressure is P_gauge = P_pipe - P_atm. Rearranging the equation:\n"
        f"P_gauge = (rho2 * g * h2) - (rho1 * g * h1)\n\n"

        f"**Step 4:** Calculate the Individual Pressure Terms\n"
        f"Pressure from manometer fluid column = {rho2} * {GRAVITY} * {h2} = {pressure_term2_pa:.2f} Pa\n"
        f"Pressure from pipe fluid column = {rho1} * {GRAVITY} * {h1} = {pressure_term1_pa:.2f} Pa\n\n"

        f"**Step 5:** Calculate the Final Gauge Pressure\n"
        f"P_gauge = {pressure_term2_pa:.2f} Pa - {pressure_term1_pa:.2f} Pa = {gauge_pressure_pa:.2f} Pa\n\n"

        f"**Step 6:** Convert the Answer to Kilopascals (kPa)\n"
        f"P_gauge = {gauge_pressure_pa:.2f} Pa / 1000 = {round(gauge_pressure_kpa, precision)} kPa\n\n"

        f"**Answer:**\n"
        f"The gauge pressure in the pipe is {round(gauge_pressure_kpa, precision)} kPa."
    )

    return question, solution


# Template 4 (Intermediate)
def template_floating_object_submersion_depth():
    """
    Fluid Statics: Floating Object Submersion Depth

    Scenario:
        This template applies the principle of buoyancy to determine how deep
        an object with a uniform cross-section (like a rectangular block or
        cylinder) will float in a liquid. The solution requires equating the
        object's total weight to the buoyant force acting on its submerged part.

    Core Equations:
        - Weight (W) = rho_object * g * V_total
        - Buoyant Force (F_B) = rho_fluid * g * V_submerged
        - At equilibrium (floating): W = F_B
        - This simplifies to: rho_object * V_total = rho_fluid * V_submerged

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the submersion depth of a floating object.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs to ensure a valid floating scenario

    # Randomly select an object material and a fluid
    obj_material, rho_object = random.choice(list(MATERIAL_DENSITIES.items()))
    fluid_name, rho_fluid = random.choice(list(FLUID_DENSITIES.items()))

    # CRITICAL: Ensure the object will actually float ---
    # Re-select until the object's density is less than the fluid's density.
    max_attempts = 20
    attempt = 0
    while rho_object >= rho_fluid and attempt < max_attempts:
        obj_material, rho_object = random.choice(list(MATERIAL_DENSITIES.items()))
        fluid_name, rho_fluid = random.choice(list(FLUID_DENSITIES.items()))
        attempt += 1
    
    # Fallback in the rare case a valid pair isn't found quickly
    if rho_object >= rho_fluid:
        obj_material, rho_object = "Pine Wood", 500
        fluid_name, rho_fluid = "Fresh Water", 998

    # Randomly choose the object's shape (uniform cross-section)
    shape = random.choice(["rectangular block", "cylinder"])

    # Randomize dimensions based on shape
    if shape == "rectangular block":
        length = round(random.uniform(0.5, 3.0), 2)
        width = round(random.uniform(0.2, 2.0), 2)
        height = round(random.uniform(0.1, 1.5), 2) # This is the total vertical height
        shape_dims_str = f"dimensions {length} m (length) x {width} m (width) x {height} m (height)"
        shape_dims_given_str = (f"  - Length (L): {length} m\n"
                                f"  - Width (W): {width} m\n"
                                f"  - Total Height (H): {height} m")
    else:  # shape == "cylinder"
        radius = round(random.uniform(0.1, 1.5), 2)
        height = round(random.uniform(0.2, 2.5), 2) # This is the total vertical height
        shape_dims_str = f"a radius of {radius} m and a total height of {height} m"
        shape_dims_given_str = (f"  - Radius (r): {radius} m\n"
                                f"  - Total Height (H): {height} m")

    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations for the solution

    # The core insight is that for a uniform cross-section, the area cancels out.
    # W = F_B  =>  rho_obj * g * (Area * H) = rho_fluid * g * (Area * h_sub)
    # Simplifying gives: rho_obj * H = rho_fluid * h_sub
    submersion_depth = (rho_object / rho_fluid) * height

    # 3. Generate the question and solution strings

    question = (
        f"A {shape} made of {obj_material.lower()} (density = {rho_object} kg/m^3) has "
        f"{shape_dims_str}. The object is placed in a large tank of {fluid_name.lower()} "
        f"(density = {rho_fluid} kg/m^3). Assuming the object floats in a stable, upright "
        f"position, calculate its submersion depth (the vertical height of the object "
        f"that is below the fluid surface)."
    )

    solution = (
        f"**Given:**\n"
        f"Object Shape: {shape.capitalize()}\n"
        f"{shape_dims_given_str}\n"
        f"Object Density (rho_obj): {rho_object} kg/m^3\n"
        f"Fluid Density (rho_fluid): {rho_fluid} kg/m^3\n\n"

        f"**Step 1:** State the Principle of Flotation\n"
        f"For an object to float, its total weight (W) must be equal to the buoyant force (F_B) "
        f"exerted by the fluid. The buoyant force is the weight of the displaced fluid.\n"
        f"Equilibrium Condition: W = F_B\n\n"

        f"**Step 2:** Express Weight and Buoyant Force using Densities\n"
        f"Weight (W) = rho_obj * g * V_total\n"
        f"Buoyant Force (F_B) = rho_fluid * g * V_submerged\n\n"
        f"Setting them equal: \n"
        f"rho_obj * g * V_total = rho_fluid * g * V_submerged\n\n"

        f"**Step 3:** Simplify for a Uniform Cross-Section\n"
        f"The acceleration of gravity (g) cancels from both sides. For an object with a uniform "
        f"cross-sectional area (A), the volumes can be expressed as V = A * height.\n"
        f"V_total = A * H_total\n"
        f"V_submerged = A * h_submerged\n\n"
        f"Substituting these into the equation:\n"
        f"rho_obj * (A * H_total) = rho_fluid * (A * h_submerged)\n\n"
        f"The cross-sectional area (A) also cancels out, leaving a simple ratio:\n"
        f"rho_obj * H_total = rho_fluid * h_submerged\n\n"

        f"**Step 4:** Solve for the Submersion Depth (h_submerged)\n"
        f"Rearranging the formula to solve for the unknown depth:\n"
        f"h_submerged = (rho_obj / rho_fluid) * H_total\n\n"
        f"Plugging in the given values:\n"
        f"h_submerged = ({rho_object} / {rho_fluid}) * {height}\n"
        f"h_submerged = {round(submersion_depth, precision)} m\n\n"

        f"**Answer:**\n"
        f"The submersion depth of the {shape} is {round(submersion_depth, precision)} m."
    )

    return question, solution


# Template 5 (Advanced)
def template_hydrostatic_force_on_plane():
    """
    Fluid Statics: Hydrostatic Force and Center of Pressure on a Submerged Plane

    Scenario:
        This template generates a comprehensive problem requiring the calculation of
        the magnitude of the resultant hydrostatic force on a submerged plane
        surface (e.g., a gate, window) and the location where this force acts,
        known as the center of pressure.

    Core Equations:
        - Resultant Force (F_R) = rho * g * h_c * A
        - Center of Pressure (y_R) = y_c + (I_xc / (y_c * A))
          where:
            - rho = fluid density
            - g = acceleration of gravity
            - h_c = vertical depth from free surface to the centroid of the area
            - A = area of the submerged surface
            - y_c = inclined distance from free surface to the centroid
            - I_xc = area moment of inertia about the centroidal axis parallel to the surface

    Returns:
        tuple: A tuple containing:
            - str: A question about hydrostatic force and center of pressure.
            - str: A detailed, step-by-step solution.
    """
    # 1. Parameterize the inputs with random values

    # Select fluid and shape
    fluid_name, rho = random.choice(list(FLUID_DENSITIES.items()))
    shape = random.choice(["rectangle", "circle"])

    # Set geometric and positional parameters
    angle_deg = random.choice([30, 45, 60, 90]) # Angle with the horizontal
    angle_rad = math.radians(angle_deg)
    h_top = round(random.uniform(0.5, 5.0), 2) # Vertical depth to top edge of the plane

    # Randomize dimensions based on shape and calculate geometric properties
    if shape == "rectangle":
        b = round(random.uniform(1.0, 4.0), 2)  # width
        h = round(random.uniform(1.0, 4.0), 2)  # height
        area = b * h
        I_xc = (b * h**3) / 12
        dist_to_centroid_along_plane = h / 2
        shape_desc = f"a rectangular gate with a width of {b} m and a height of {h} m"
        shape_given = f"  - Shape: Rectangle (width b = {b} m, height h = {h} m)"
    else:  # shape == "circle"
        r = round(random.uniform(0.5, 2.0), 2) # radius
        area = math.pi * r**2
        I_xc = (math.pi * r**4) / 4
        dist_to_centroid_along_plane = r
        shape_desc = f"a circular viewport with a radius of {r} m"
        shape_given = f"  - Shape: Circle (radius r = {r} m)"

    precision = 4

    # 2. Perform the core calculations for the solution

    # Step A: Determine the position of the centroid (yc and hc)
    # yc is the inclined distance from the free surface to the centroid
    # hc is the vertical depth from the free surface to the centroid
    if angle_deg == 90:
        y_top = h_top
        angle_desc = "is positioned vertically"
    else:
        y_top = h_top / math.sin(angle_rad)
        angle_desc = f"is inclined at an angle of {angle_deg} degrees to the horizontal"

    y_c = y_top + dist_to_centroid_along_plane
    h_c = y_c * math.sin(angle_rad)

    # Step B: Calculate the resultant force (F_R)
    force_resultant_N = rho * GRAVITY * h_c * area
    force_resultant_kN = force_resultant_N / 1000

    # Step C: Calculate the location of the center of pressure (y_R)
    # y_R is the inclined distance from the free surface to the center of pressure
    y_R = y_c + (I_xc / (y_c * area))

    # 3. Generate the question and solution strings

    question = (
        f"A submerged {shape_desc} {angle_desc} in a tank of {fluid_name.lower()} "
        f"(density = {rho} kg/m^3). The top edge of the gate is {h_top} m vertically below "
        f"the free surface. Calculate:\n"
        f"  a) The magnitude of the resultant hydrostatic force on the gate.\n"
        f"  b) The location of the center of pressure, measured along the incline of the gate "
        f"from the free surface."
    )

    solution = (
        f"**Given:**\n"
        f"{shape_given}\n"
        f"Fluid: {fluid_name} (rho = {rho} kg/m^3)\n"
        f"Vertical depth to top edge (h_top): {h_top} m\n"
        f"Angle of inclination (theta): {angle_deg} degrees\n\n"

        f"**Step 1:** Calculate Geometric Properties of the Gate\n"
        f"Area (A): {area:.{precision}f} m^2\n"
        f"Moment of Inertia about centroid (I_xc): {I_xc:.{precision}f} m^4\n\n"

        f"**Step 2:** Determine the Location of the Centroid (yc and hc)\n"
        f"The centroid is the geometric center of the gate. We need its position relative to the free surface.\n"
        f"  - 'y' distances are measured along the plane's incline.\n"
        f"  - 'h' distances are measured vertically.\n\n"
        f"Inclined distance from surface to top edge (y_top) = h_top / sin(theta)\n"
        f"  y_top = {h_top} / sin({angle_deg}°) = {y_top:.{precision}f} m\n"
        f"Distance from top edge to centroid along plane = {dist_to_centroid_along_plane:.{precision}f} m\n"
        f"Inclined distance from surface to centroid (y_c) = y_top + dist_to_centroid\n"
        f"  y_c = {y_top:.{precision}f} + {dist_to_centroid_along_plane:.{precision}f} = {y_c:.{precision}f} m\n\n"
        f"Vertical depth to centroid (h_c) = y_c * sin(theta)\n"
        f"  h_c = {y_c:.{precision}f} * sin({angle_deg}°) = {h_c:.{precision}f} m\n\n"

        f"**Step 3:** Calculate the Resultant Hydrostatic Force (F_R)\n"
        f"The force is the pressure at the centroid multiplied by the total area.\n"
        f"F_R = rho * g * h_c * A\n"
        f"F_R = {rho} * {GRAVITY} * {h_c:.{precision}f} * {area:.{precision}f}\n"
        f"F_R = {force_resultant_N:.2f} N\n"
        f"F_R = {force_resultant_kN:.{precision}f} kN\n\n"

        f"**Step 4:** Calculate the Center of Pressure (y_R)\n"
        f"The center of pressure is the point where the resultant force acts. It is always below the centroid.\n"
        f"y_R = y_c + (I_xc / (y_c * A))\n"
        f"y_R = {y_c:.{precision}f} + ({I_xc:.{precision}f} / ({y_c:.{precision}f} * {area:.{precision}f}))\n"
        f"y_R = {y_R:.{precision}f} m\n\n"

        f"**Answer:**\n"
        f"a) The magnitude of the resultant hydrostatic force is **{force_resultant_kN:.{precision}f} kN**.\n"
        f"b) The center of pressure is located **{y_R:.{precision}f} m** from the free surface, measured down along the angle of the gate."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each fluid statics template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/fluid_mechanics/fluid_statics.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_hydrostatic_pressure_at_depth, "hydrostatic_pressure_at_depth", "Easy"),
        (template_basic_buoyant_force, "basic_buoyant_force", "Easy"),
        (template_utube_manometer, "utube_manometer", "Intermediate"),
        (template_floating_object_submersion_depth, "floating_object_submersion_depth", "Intermediate"),
        (template_hydrostatic_force_on_plane, "hydrostatic_force_on_plane", "Advanced"),
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
