import random
import math
from data.templates.branches.mechanical_engineering.constants import SHEAR_MODULUS_VALUES

# Template 1 (Easy)
def template_shear_stress_torsion():
    """
    Torsion: Shear Stress and Polar Moment of Inertia

    Scenario:
        This template generates a foundational problem testing the ability to calculate
        the polar moment of inertia (J) for a solid or hollow circular shaft.
        It then uses this value to determine the maximum shearing stress (tau_max)
        resulting from an applied torque.

    Core Equations:
        tau_max = (T * c) / J
        J_solid = (pi / 2) * c^4
        J_hollow = (pi / 2) * (c_outer^4 - c_inner^4)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the maximum shearing stress in a shaft.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    torque = round(random.uniform(100.0, 5000.0), 2)  # Torque in N.m
    shaft_type = random.choice(['solid', 'hollow'])
    
    # Ensure diameters are distinct and practical
    d_outer = random.randint(30, 150)  # Outer diameter in mm
    
    # For hollow shafts, ensure the inner diameter is smaller than the outer
    if shaft_type == 'hollow':
        # Ensure inner diameter is at least 10mm smaller to be meaningful
        max_inner_dia = max(20, d_outer - 10) 
        d_inner = random.randint(20, max_inner_dia)
    else:
        d_inner = 0

    # Standardize precision for all calculations and outputs
    precision = 3

    # 2. Perform the core calculations for the solution
    
    # Step A: Convert diameters to radii in meters
    c_outer = d_outer / 2000.0  # Convert mm to m and get radius
    if shaft_type == 'hollow':
        c_inner = d_inner / 2000.0
    
    # Step B: Calculate the polar moment of inertia (J)
    if shaft_type == 'solid':
        polar_moment_J = (math.pi / 2) * (c_outer ** 4)
        j_calculation_str = f"J = (pi / 2) * c^4 = (pi / 2) * ({c_outer})^4 = {polar_moment_J:.3e} m^4"
    else: # shaft_type == 'hollow'
        polar_moment_J = (math.pi / 2) * (c_outer ** 4 - c_inner ** 4)
        j_calculation_str = (
            f"J = (pi / 2) * (c_outer^4 - c_inner^4)\n"
            f"   J = (pi / 2) * (({c_outer})^4 - ({c_inner})^4) = {polar_moment_J:.3e} m^4"
        )
        
    # Step C: Apply the torsion formula to find the stress in Pascals
    # Note: 'c' in the formula refers to the outermost radius, c_outer.
    tau_max_pascals = (torque * c_outer) / polar_moment_J
    
    # Step D: Convert the final answer to megapascals (MPa)
    tau_max_mpa = tau_max_pascals / 1e6

    # 3. Generate the question and solution strings
    
    # Construct the part of the question describing the geometry
    if shaft_type == 'solid':
        geometry_desc = f"a solid circular shaft with an outer diameter of {d_outer} mm"
    else: # shaft_type == 'hollow'
        geometry_desc = (
            f"a hollow circular shaft with an outer diameter of {d_outer} mm "
            f"and an inner diameter of {d_inner} mm"
        )

    question = (
        f"A {geometry_desc} is subjected to a torque of {torque} N.m. "
        f"Determine the maximum shearing stress in the shaft."
    )

    solution = (
        f"**Given:**\n"
        f"Torque (T) = {torque} N.m\n"
        f"Shaft Type = {shaft_type.capitalize()}\n"
        f"Outer Diameter (d_outer) = {d_outer} mm\n"
    )
    if shaft_type == 'hollow':
        solution += f"  - Inner Diameter (d_inner) = {d_inner} mm\n\n"
    else:
        solution += "\n"

    solution += (
        f"**Step 1:** Convert diameters to radii and express in meters.\n"
        f"The maximum stress occurs at the outer surface, so we use the outer radius (c) for the stress calculation.\n"
        f"Outer radius (c) = d_outer / 2 = {d_outer} / 2 = {d_outer/2.0} mm = {c_outer} m\n"
    )
    if shaft_type == 'hollow':
        solution += f"  - Inner radius (c_inner) = d_inner / 2 = {d_inner} / 2 = {d_inner/2.0} mm = {c_inner} m\n"
    solution += "\n"
    
    solution += (
        f"**Step 2:** Calculate the polar moment of inertia (J) for the shaft's cross-section.\n"
        f"For a {shaft_type} shaft:\n"
        f"  {j_calculation_str}\n\n"
        
        f"**Step 3:** Apply the torsion formula to find the maximum shearing stress (tau_max).\n"
        f"The formula is: tau_max = (T * c) / J\n"
        f"T = {torque} N.m\n"
        f"c = {c_outer} m\n"
        f"J = {polar_moment_J:.3e} m^4\n"
        f"tau_max = ({torque} * {c_outer}) / {polar_moment_J:.3e}\n"
        f"tau_max = {tau_max_pascals:.3e} Pa\n\n"

        f"**Step 4:** Convert the stress from Pascals (Pa) to Megapascals (MPa).\n"
        f"1 MPa = 1,000,000 Pa\n"
        f"tau_max = {tau_max_pascals:.3e} Pa / 1e6 = {round(tau_max_mpa, precision)} MPa\n\n"

        f"**Answer:**\n"
        f"The maximum shearing stress in the shaft is {round(tau_max_mpa, precision)} MPa."
    )

    return question, solution


# Template 2 (Easy)
def template_angle_of_twist():
    """
    Torsion: Angle of Twist Calculation

    Scenario:
        This template assesses the ability to calculate the total angle of twist (phi)
        for a uniform solid circular shaft. It requires a correct understanding of the
        relationship between torque, length, material properties (Shear Modulus), and
        the shaft's geometry (Polar Moment of Inertia).

    Core Equations:
        phi = (T * L) / (J * G)
        J_solid = (pi / 2) * c^4

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the angle of twist in a shaft.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    torque = round(random.uniform(500.0, 6000.0), 1)   # Torque in N.m
    length = round(random.uniform(0.5, 4.0), 2)       # Length in m
    diameter = random.randint(25, 100)                # Diameter in mm
    
    # Randomly select a material and its properties
    material_name, shear_modulus_gpa = random.choice(list(SHEAR_MODULUS_VALUES.items()))
    
    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations for the solution
    
    # Step A: Convert diameter (mm) to radius (m)
    c_radius_m = diameter / 2000.0
    
    # Step B: Calculate the polar moment of inertia (J)
    polar_moment_J = (math.pi / 2) * (c_radius_m ** 4)
    
    # Step C: Ensure all units are consistent (convert G from GPa to Pa)
    shear_modulus_pa = shear_modulus_gpa * 1e9
    
    # Step D: Apply the angle of twist formula to find the angle in radians
    angle_rad = (torque * length) / (polar_moment_J * shear_modulus_pa)
    
    # Step E: Convert the result from radians to degrees
    angle_deg = math.degrees(angle_rad)

    # 3. Generate the question and solution strings
    
    question = (
        f"A solid {material_name.lower()} shaft with a diameter of {diameter} mm and a length of {length} m "
        f"is subjected to a torque of {torque} N.m. "
        f"Given that the shear modulus (G) for {material_name.lower()} is {shear_modulus_gpa} GPa, "
        f"calculate the total angle of twist. Provide the answer in both radians and degrees."
    )

    solution = (
        f"**Given:**\n"
        f"Torque (T) = {torque} N.m\n"
        f"Length (L) = {length} m\n"
        f"Diameter (d) = {diameter} mm\n"
        f"Material = {material_name}\n"
        f"Shear Modulus (G) = {shear_modulus_gpa} GPa\n\n"

        f"**Step 1:** Convert the diameter to a radius in meters.\n"
        f"Radius (c) = d / 2 = {diameter} / 2 = {diameter/2.0} mm = {c_radius_m} m\n\n"

        f"**Step 2:** Calculate the polar moment of inertia (J) for the solid circular shaft.\n"
        f"Formula: J = (pi / 2) * c^4\n"
        f"J = (pi / 2) * ({c_radius_m})^4 = {polar_moment_J:.4e} m^4\n\n"
        
        f"**Step 3:** Ensure consistent units for the angle of twist calculation.\n"
        f"The shear modulus must be in Pascals (Pa) to be consistent with N and m.\n"
        f"G = {shear_modulus_gpa} GPa = {shear_modulus_pa:.2e} Pa\n\n"

        f"**Step 4:** Apply the angle of twist formula to find the angle in radians.\n"
        f"Formula: phi = (T * L) / (J * G)\n"
        f"phi = ({torque} * {length}) / ({polar_moment_J:.4e} * {shear_modulus_pa:.2e})\n"
        f"phi = {round(angle_rad, precision)} radians\n\n"

        f"**Step 5:** Convert the angle from radians to degrees.\n"
        f"Angle in degrees = Angle in radians * (180 / pi)\n"
        f"Angle = {round(angle_rad, precision)} * (180 / pi) = {round(angle_deg, precision)} degrees\n\n"

        f"**Answer:**\n"
        f"The total angle of twist is {round(angle_rad, precision)} radians, which is equivalent to "
        f"{round(angle_deg, precision)} degrees."
    )

    return question, solution


# Template 3 (Intermediate)
def template_shaft_design_power():
    """
    Torsion: Power Transmission and Shaft Design

    Scenario:
        This template generates a classic design problem. It tests the ability to
        first determine the torque on a shaft from its power and rotational speed,
        and then use that torque to calculate the minimum required shaft diameter
        based on an allowable shearing stress.

    Core Equations:
        P = 2 * pi * f * T
        tau_allow = (T * c) / J
        For a solid shaft, this simplifies to: c = ( (2 * T) / (pi * tau) )^(1/3)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the minimum shaft diameter.
            - str: A step-by-step solution to the design problem.
    """
    # 1. Parameterize the inputs with random values
    power_kw = round(random.uniform(10.0, 500.0), 1)
    allowable_stress_mpa = random.randint(40, 120)
    
    # Randomly choose to provide frequency in Hz or RPM to add variety
    use_rpm = random.choice([True, False])
    if use_rpm:
        frequency_rpm = random.randint(500, 3000)
        frequency_hz = frequency_rpm / 60.0
        frequency_str = f"{frequency_rpm} RPM"
    else:
        frequency_hz = random.randint(10, 100)
        frequency_rpm = frequency_hz * 60
        frequency_str = f"{frequency_hz} Hz"

    # Standardize precision for final outputs
    precision = 2

    # 2. Perform the core calculations for the solution
    
    # Step A: Convert primary units to base SI units (W, Pa, Hz)
    power_w = power_kw * 1000
    allowable_stress_pa = allowable_stress_mpa * 1e6
    
    # Step B: Calculate the torque (T) on the shaft using the power formula
    # P = 2 * pi * f * T  =>  T = P / (2 * pi * f)
    torque = power_w / (2 * math.pi * frequency_hz)
    
    # Step C: Calculate the required radius (c) using the rearranged torsion formula
    # tau = (T*c) / J = (T*c) / (pi/2 * c^4) = 2*T / (pi*c^3)
    # c^3 = (2 * T) / (pi * tau)
    c_cubed = (2 * torque) / (math.pi * allowable_stress_pa)
    c_radius_m = c_cubed ** (1/3)
    
    # Step D: Calculate the diameter in meters and then convert to millimeters
    d_diameter_m = c_radius_m * 2
    d_diameter_mm = d_diameter_m * 1000

    # 3. Generate the question and solution strings
    
    question = (
        f"A motor is required to transmit {power_kw} kW of power at a rotational speed of {frequency_str}. "
        f"If the solid circular shaft is to be made from a material with an allowable shearing stress of {allowable_stress_mpa} MPa, "
        f"determine the minimum required diameter for the shaft."
    )

    solution = (
        f"**Given:**\n"
        f"Power (P) = {power_kw} kW\n"
        f"Rotational Speed = {frequency_str}\n"
        f"Allowable Shearing Stress (tau_allow) = {allowable_stress_mpa} MPa\n\n"

        f"**Step 1:** Convert the given values to base SI units (Watts, Pascals, Hertz).\n"
        f"Power (P) = {power_kw} kW = {power_w} W\n"
        f"Allowable Stress (tau_allow) = {allowable_stress_mpa} MPa = {allowable_stress_pa:.1e} Pa\n"
    )

    if use_rpm:
        solution += (
            f"  - Frequency (f) = {frequency_rpm} RPM = {frequency_rpm} / 60 = {round(frequency_hz, 2)} Hz\n\n"
        )
    else:
        solution += f"  - Frequency (f) = {frequency_hz} Hz\n\n"

    solution += (
        f"**Step 2:** Calculate the torque (T) exerted on the shaft.\n"
        f"The relationship between power, torque, and frequency is P = 2 * pi * f * T.\n"
        f"Rearranging for torque: T = P / (2 * pi * f)\n"
        f"T = {power_w} / (2 * pi * {round(frequency_hz, 2)}) = {round(torque, 2)} N.m\n\n"
        
        f"**Step 3:** Determine the required shaft radius (c) using the torsion formula.\n"
        f"The formula for maximum stress in a solid shaft is tau = (T * c) / J, where J = (pi/2) * c^4.\n"
        f"This simplifies to tau = 2 * T / (pi * c^3).\n"
        f"Rearranging to solve for the radius: c^3 = (2 * T) / (pi * tau_allow)\n"
        f"c^3 = (2 * {round(torque, 2)}) / (pi * {allowable_stress_pa:.1e}) = {c_cubed:.3e} m^3\n"
        f"c = ({c_cubed:.3e})^(1/3) = {round(c_radius_m, 4)} m\n\n"

        f"**Step 4:** Calculate the minimum diameter from the radius.\n"
        f"Diameter (d) = 2 * c\n"
        f"d = 2 * {round(c_radius_m, 4)} = {round(d_diameter_m, 4)} m\n"
        f"In millimeters, d = {round(d_diameter_mm, precision)} mm\n\n"
        
        f"**Answer:**\n"
        f"The minimum required diameter for the solid circular shaft is {round(d_diameter_mm, precision)} mm."
    )

    return question, solution


# Template 4 (Intermediate)
def template_composite_shafts_series():
    """
    Torsion: Composite Shafts in Series

    Scenario:
        This problem involves a shaft made of two different segments joined end-to-end.
        It tests the understanding that the total angle of twist at the free end is
        the algebraic sum of the angles of twist of each individual segment.

    Core Equations:
        phi_total = sum( (T_i * L_i) / (J_i * G_i) )
        For this case: phi_total = phi_1 + phi_2
        J_solid = (pi / 2) * c^4

    Returns:
        tuple: A tuple containing:
            - str: A question about a composite shaft in series.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    torque = round(random.uniform(200.0, 7500.0), 1)

    # Properties for Segment 1 (AB)
    l1 = round(random.uniform(0.5, 2.5), 2)
    d1 = random.randint(40, 120)
    mat1_name, g1_gpa = random.choice(list(SHEAR_MODULUS_VALUES.items()))

    # Properties for Segment 2 (BC)
    l2 = round(random.uniform(0.5, 2.5), 2)
    d2 = random.randint(30, d1) # Ensure d2 is not larger than d1 for a typical stepped shaft
    mat2_name, g2_gpa = random.choice(list(SHEAR_MODULUS_VALUES.items()))
    
    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations
    
    # --- Calculations for Segment 1 (AB) ---
    c1_m = d1 / 2000.0
    g1_pa = g1_gpa * 1e9
    j1 = (math.pi / 2) * (c1_m ** 4)
    phi1_rad = (torque * l1) / (j1 * g1_pa)

    # --- Calculations for Segment 2 (BC) ---
    c2_m = d2 / 2000.0
    g2_pa = g2_gpa * 1e9
    j2 = (math.pi / 2) * (c2_m ** 4)
    phi2_rad = (torque * l2) / (j2 * g2_pa)
    
    # --- Total Angle of Twist ---
    phi_total_rad = phi1_rad + phi2_rad
    phi_total_deg = math.degrees(phi_total_rad)

    # 3. Generate the question and solution strings
    
    question = (
        f"A composite shaft consists of two segments, AB and BC, rigidly connected at B. "
        f"Segment AB is a solid {mat1_name.lower()} shaft of length {l1} m and diameter {d1} mm. "
        f"Segment BC is a solid {mat2_name.lower()} shaft of length {l2} m and diameter {d2} mm. "
        f"The shaft is fixed at end A, and a torque of {torque} N.m is applied at the free end C. "
        f"Calculate the total angle of twist at end C. "
        f"(Use G = {g1_gpa} GPa for {mat1_name.lower()} and G = {g2_gpa} GPa for {mat2_name.lower()})."
    )

    solution = (
        f"**Given:**\n"
        f"Applied Torque (T) = {torque} N.m\n"
        f"Segment AB: L1={l1} m, d1={d1} mm, Material={mat1_name} (G1={g1_gpa} GPa)\n"
        f"Segment BC: L2={l2} m, d2={d2} mm, Material={mat2_name} (G2={g2_gpa} GPa)\n\n"
        
        f"**Step 1:** Analyze the System\n"
        f"The shaft is fixed at A, so the torque T = {torque} N.m is transmitted through both segments AB and BC. "
        f"The total angle of twist at C is the sum of the twist in segment AB and the twist in segment BC.\n"
        f"phi_total = phi_AB + phi_BC\n\n"

        f"**Step 2:** Calculate Angle of Twist for Segment AB (phi_AB)\n"
        f"Convert units: c1 = {d1}/2 mm = {c1_m} m; G1 = {g1_gpa} GPa = {g1_pa:.2e} Pa\n"
        f"Polar Moment of Inertia (J1) = (pi/2) * c1^4 = (pi/2) * ({c1_m})^4 = {j1:.3e} m^4\n"
        f"Angle of Twist (phi_AB) = (T * L1) / (J1 * G1)\n"
        f"phi_AB = ({torque} * {l1}) / ({j1:.3e} * {g1_pa:.2e}) = {round(phi1_rad, precision + 1)} radians\n\n"
        
        f"**Step 3:** Calculate Angle of Twist for Segment BC (phi_BC)\n"
        f"Convert units: c2 = {d2}/2 mm = {c2_m} m; G2 = {g2_gpa} GPa = {g2_pa:.2e} Pa\n"
        f"Polar Moment of Inertia (J2) = (pi/2) * c2^4 = (pi/2) * ({c2_m})^4 = {j2:.3e} m^4\n"
        f"Angle of Twist (phi_BC) = (T * L2) / (J2 * G2)\n"
        f"phi_BC = ({torque} * {l2}) / ({j2:.3e} * {g2_pa:.2e}) = {round(phi2_rad, precision + 1)} radians\n\n"
        
        f"**Step 4:** Calculate Total Angle of Twist at End C\n"
        f"phi_total = phi_AB + phi_BC\n"
        f"phi_total = {round(phi1_rad, precision + 1)} + {round(phi2_rad, precision + 1)} = {round(phi_total_rad, precision)} radians\n"
        f"To convert to degrees: Angle_deg = Angle_rad * (180 / pi)\n"
        f"phi_total = {round(phi_total_rad, precision)} * (180 / pi) = {round(phi_total_deg, precision)} degrees\n\n"
        
        f"**Answer:**\n"
        f"The total angle of twist at the free end C is {round(phi_total_rad, precision)} radians, "
        f"or {round(phi_total_deg, precision)} degrees."
    )

    return question, solution


# Template 5 (Advanced)
def template_statically_indeterminate_shaft():
    """
    Torsion: Statically Indeterminate Shaft

    Scenario:
        This problem deals with a shaft that is fixed at both ends and has a torque
        applied at an intermediate point. Since there are two unknown reaction torques
        but only one static equilibrium equation, the problem is statically
        indeterminate. It must be solved by considering both static equilibrium and
        the geometric compatibility of deformation (i.e., the total angle of twist is zero).

    Core Equations:
        1. Statics: T_A + T_B = T_applied
        2. Compatibility: phi_AC + phi_CB = 0  =>  (T_A * L_AC) / (JG) = (T_B * L_BC) / (JG)
           This simplifies to: T_A * L_AC = T_B * L_BC

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the reaction torques.
            - str: A step-by-step solution using statics and compatibility.
    """
    # 1. Parameterize the inputs with random values
    total_length = round(random.uniform(2.0, 6.0), 2)
    applied_torque = round(random.uniform(1000.0, 15000.0), -2)
    diameter = random.randint(50, 150)
    
    # Ensure the torque is applied at a non-trivial location
    pos_L_AC = round(random.uniform(0.2 * total_length, 0.8 * total_length), 2)
    pos_L_BC = total_length - pos_L_AC

    material_name, shear_modulus_gpa = random.choice(list(SHEAR_MODULUS_VALUES.items()))
    
    # Standardize precision for final outputs
    precision = 2

    # 2. Perform the core calculations
    # Based on the system of equations:
    # (1) T_A + T_B = T_applied
    # (2) T_A * L_AC = T_B * L_BC  =>  T_B = T_A * (L_AC / L_BC)
    # Substitute (2) into (1):
    # T_A + T_A * (L_AC / L_BC) = T_applied
    # T_A * (1 + L_AC/L_BC) = T_applied
    # T_A * ( (L_BC + L_AC) / L_BC ) = T_applied
    # T_A * (L / L_BC) = T_applied
    # T_A = T_applied * (L_BC / L)
    
    reaction_torque_A = applied_torque * (pos_L_BC / total_length)
    reaction_torque_B = applied_torque * (pos_L_AC / total_length)

    # 3. Generate the question and solution strings
    
    question = (
        f"A solid {material_name.lower()} shaft of diameter {diameter} mm and length {total_length} m is "
        f"fixed at both ends, A and B. A torque of {int(applied_torque)} N.m is applied at point C, "
        f"located {pos_L_AC} m from end A. The shear modulus for the material is {shear_modulus_gpa} GPa. "
        f"Determine the reaction torques at the fixed supports, T_A and T_B."
    )

    solution = (
        f"**Given:**\n"
        f"Total Length (L) = {total_length} m\n"
        f"Applied Torque (T_applied) = {int(applied_torque)} N.m\n"
        f"Location of Torque from A (L_AC) = {pos_L_AC} m\n"
        f"Diameter (d) = {diameter} mm, Material = {material_name}\n\n"
        
        f"**Analysis:**\n"
        f"This problem is statically indeterminate because there are two unknown reaction torques (T_A and T_B) "
        f"and only one relevant equation from statics. We need an additional equation from the deformation of the shaft.\n\n"

        f"**Step 1:** Statics Equilibrium Equation\n"
        f"For the shaft to be in rotational equilibrium, the sum of all torques must be zero. Let's assume T_A and T_B act in the opposite direction to T_applied.\n"
        f"(1) T_A + T_B = T_applied = {int(applied_torque)} N.m\n\n"

        f"**Step 2:** Compatibility Equation\n"
        f"Since the shaft is fixed at both ends, the total angle of twist from A to B must be zero. The twist from A to C and the twist from C to B must cancel each other out.\n"
        f"phi_A_to_B = phi_AC + phi_CB = 0\n"
        f"The torque in section AC is T_A. The torque in section CB is T_applied - T_A = T_B. For our compatibility equation, it is easier to think of the torque in CB as T_B acting from the other direction.\n"
        f"phi_AC = (T_A * L_AC) / (J*G)\n"
        f"phi_CB = -(T_B * L_BC) / (J*G)  (The negative sign indicates it twists in the opposite direction)\n"
        f"(T_A * L_AC) / (J*G) - (T_B * L_BC) / (J*G) = 0\n"
        f"The terms J and G are constant for the shaft and cancel out, leaving a relationship between the torques and lengths:\n"
        f"(2) T_A * L_AC = T_B * L_BC\n\n"

        f"**Step 3:** Solve the System of Two Equations\n"
        f"We have two equations:\n"
        f"(1) T_A + T_B = {int(applied_torque)}\n"
        f"(2) T_A * {pos_L_AC} = T_B * {round(pos_L_BC, 2)}\n"
        f"From equation (2), we can express T_B in terms of T_A:\n"
        f"T_B = T_A * ({pos_L_AC} / {round(pos_L_BC, 2)})\n"
        f"Substitute this into equation (1):\n"
        f"T_A + T_A * ({pos_L_AC} / {round(pos_L_BC, 2)}) = {int(applied_torque)}\n"
        f"T_A * (1 + {round(pos_L_AC / pos_L_BC, 3)}) = {int(applied_torque)}\n"
        f"T_A = {int(applied_torque)} / {round(1 + (pos_L_AC / pos_L_BC), 3)} = {round(reaction_torque_A, precision)} N.m\n"
        f"Now find T_B using equation (1):\n"
        f"T_B = {int(applied_torque)} - T_A = {int(applied_torque)} - {round(reaction_torque_A, precision)} = {round(reaction_torque_B, precision)} N.m\n\n"

        f"**Answer:**\n"
        f"The reaction torques at the supports are:\n"
        f"T_A = {round(reaction_torque_A, precision)} N.m\n"
        f"T_B = {round(reaction_torque_B, precision)} N.m"
    )

    return question, solution


def main():
    """
    Generate numerous instances of each torsion template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/mechanics_of_materials/torsion.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_shear_stress_torsion, "shear_stress_torsion", "Easy"),
        (template_angle_of_twist, "angle_of_twist", "Easy"),
        (template_shaft_design_power, "shaft_design_power", "Intermediate"),
        (template_composite_shafts_series, "composite_shafts_series", "Intermediate"),
        (template_statically_indeterminate_shaft, "statically_indeterminate_shaft", "Advanced"),
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
                "branch": "mechanical_engineering",
                "domain": "mechanics_of_materials",
                "area": "torsion",
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
